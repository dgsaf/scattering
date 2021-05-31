!>
module basis

#include "debug.h"

  ! debug compilation
  ! - <DEBUG_BASIS>: verbosity of debug statements, local to this file.
#define DEBUG_BASIS 2

  ! local debug verbosity cannot be lower than the global debug verbosity
#if ((defined DEBUG) && (defined DEBUG_BASIS) && (DEBUG > DEBUG_BASIS))
  DEBUG_BASIS = DEBUG
#endif

  use parameters

  implicit none

  ! private
  ! public setup_states, setup_radial, valid_states, valid_radial, &
  !     partial_waves, overlap, kinetic

  ! t_basis
  !
  ! Laguerre basis, consisting of basis states
  ! > {|phi_{i}>} for i = 1, .., <n_basis>
  ! with coordinate-space representation
  ! > phi_{i}(r, theta, phi) = (varphi_{k_{i}, l_{i}}(r) / r)
  ! >                          * Y_{l_{i}, m_{i}}(theta, phi)
  ! where
  ! > varphi_{k, l}(r) = sqrt(alpha * (k - 1)! / (k + l) * (k + 2*l)!)
  ! >                    * (2*alpha*r)^{l+1}
  ! >                    * exp(-alpha*r)
  ! >                    * L_{k - 1}^{2*l + 1}(2*alpha*r)
  ! where L_{i}^{j} are the generalised Laguerre polynomials.
  !
  ! We construct the basis for given:
  ! - <l_max>: maximum angular quantum number, <l>, considered in the basis;
  ! - <n_basis_l>: number of basis functions per <l> (from 0 to <l_max>);
  ! - <alpha_l>: value of <alpha> per <l> (from 0 to <l_max>).
  !
  ! Optionally, we may impose the following symmetries:
  ! - <m>: the magnetic quantum number;
  ! - <parity>: the parity quantum number, <parity> = (-1)^<l>;
  ! - <l>: the angular quantum number.
  !
  ! As a result of this construction, we can determine:
  ! - <n_basis>: the total number of basis states;
  ! - <k_list>: value of <k> for each basis state, k_{i} = k_list(i);
  ! - <l_list>: value of <l> for each basis state, l_{i} = l_list(i);
  ! - <m_list>: value of <m> for each basis state, m_{i} = m_list(i).
  !
  ! Furthermore, given a radial grid, <r_grid>, of length <n_r>, we also store:
  ! - <r_grid>: the radial grid;
  ! - <n_r>: the number of points in the radial grid;
  ! - <radial>: the radial basis functions, varphi_{k_{i}, l_{i}}(r), calculated
  !   on the radial grid points.
  !
  ! For the basis to be valid, we must have that:
  ! - <n_basis_l(:)> >= 0, as we cannot have a negative amount of basis states
  !   for any value of <l>;
  ! - <alpha_l(:)> > 0.0, since alpha relates to radial variable, and must be
  !   strictly positive;
  ! - <n_basis> >= 1, since we require at least one basis state in a basis;
  ! - if <m> symmetry imposed:
  !   - <l_max> >= |<m>|, as we must have <m> in {-<l>, .., <l>} for all basis
  !     states;
  ! - if <parity> symmetry imposed:
  !   - <parity> = -1, or +1;
  ! - if <l> symmetry imposed:
  !   - <l> >= 0
  !   - <l> >= |<m>|.
  !
  ! We note that the basis is ordered in the following way
  ! > {|k,l,m> : k = 1, .., n_basis(l) : l = 0, .., l_max : m = -l, .., l}
  ! where the left-most iterations are iterated through first, and where
  ! > |k,l,m> = |phi_{k_{i}, l_{i}, m_{i}}> = |phi_{i}>.
  type t_basis
    logical :: has_sym_m
    integer :: m

    logical :: has_sym_parity
    integer :: parity

    logical :: has_sym_l
    integer :: l

    integer :: l_max
    integer , allocatable :: n_basis_l(:)
    double precision , allocatable :: alpha_l(:)

    integer :: n_basis
    integer , allocatable :: k_list(:)
    integer , allocatable :: l_list(:)
    integer , allocatable :: m_list(:)

    integer :: n_r
    double precision , allocatable :: r_grid(:)
    double precision , allocatable :: radial(:, :)
  end type t_basis

contains

  ! setup_states
  !
  ! For given <basis>, <l_max>, <n_basis_l>, <alpha_l> (and optional symmetries:
  ! <m>, <parity>, <l>), the Laguerre basis is constructed (if valid).
  !
  ! The following will be calculated:
  ! - <n_basis>;
  ! - <k_list>;
  ! - <l_list>;
  ! - <m_list>.
  !
  ! Returns an error code, <i_err>, where:
  ! - 0: indicates successful execution;
  ! - 1: indicates invalid arguments.
  subroutine setup_states (basis, l_max, n_basis_l, alpha_l, i_err, m, parity, l)
    type(t_basis) , intent(inout) :: basis
    integer , intent(in) :: l_max
    integer , intent(in) :: n_basis_l(0:l_max)
    double precision , intent(in) :: alpha_l(0:l_max)
    integer , intent(out) :: i_err
    integer , optional , intent(in) :: m
    integer , optional , intent(in) :: parity
    integer , optional , intent(in) :: l
    integer :: ii, kk, ll, mm

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "subroutine setup_states()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (l_max < 0) then
      i_err = 1
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<l_max> < 0"
#endif
    else
      if (any(n_basis_l(:) < 0)) then
        i_err = 1
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "any(<n_basis_l(:)> < 0)"
#endif
      end if

      if (any(alpha_l(:) < TOL)) then
        i_err = 1
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "any(<alpha_l(:)> < TOL)"
#endif
      end if
    end if

    if (present(m)) then
      if (l_max < abs(m)) then
        i_err = 1
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "<l_max> < abs(<m>)"
#endif
      end if
    end if

    if (present(parity)) then
      if (abs(parity) /= 1) then
        i_err = 1
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "abs(<parity>) /= 1"
#endif
      end if
    end if

    if (present(l)) then
      if (l < 0) then
        i_err = 1
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "<l> < 0"
#endif
      end if

      if (l > l_max) then
        i_err = 1
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "<l> > <l_max>"
#endif
      end if

      if (present(m)) then
        if (l < abs(m)) then
          i_err = 1
#if (DEBUG_BASIS >= 2)
          write (STDERR, *) DBG, ERR, "<l> < abs(<m>)"
#endif
        end if
      end if
    end if

    ! handle invalid arguments
    if (i_err /= 0) then
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_BASIS >= 1)
      write (STDERR, *) DBG, ERR, "arguments are invalid"
#endif

      ! set scalar variables to erroneous values
      basis%has_sym_m = .false.
      basis%m = 0
      basis%has_sym_parity = .false.
      basis%parity = 0
      basis%has_sym_l = .false.
      basis%l = 0
      basis%l_max = -1
      basis%n_basis = 0
      basis%n_r = 0

#if (DEBUG_BASIS >= 1)
      write (STDERR, *) DBG, ERR, &
          "<basis> array variables will be left un-allocated"
      write (STDERR, *) DBG, ERR, &
          "<basis> scalar variables set to erroneous values"
#endif
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<basis%has_sym_m> = ", basis%has_sym_m
      write (STDERR, *) DBG, ERR, "<basis%m> = ", basis%m
      write (STDERR, *) DBG, ERR, &
          "<basis%has_sym_parity> = ", basis%has_sym_parity
      write (STDERR, *) DBG, ERR, "<basis%parity> = ", basis%parity
      write (STDERR, *) DBG, ERR, "<basis%has_sym_l> = ", basis%has_sym_l
      write (STDERR, *) DBG, ERR, "<basis%l> = ", basis%l
      write (STDERR, *) DBG, ERR, "<basis%l_max> = ", basis%l_max
      write (STDERR, *) DBG, ERR, "<basis%n_basis> = ", basis%n_basis
      write (STDERR, *) DBG, ERR, "<basis%n_r> = ", basis%n_r
#endif
#if (DEBUG_BASIS >= 1)
      write (STDERR, *) DBG, ERR, "exiting subroutine setup_states()"
#endif
      return
    end if

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "arguments are valid"
#endif

    ! if <m> symmetry imposed, set <basis%has_sym_m>, <basis%m>
    ! else, set <basis%has_sym_m>, and arbitrarily set <basis%m>
    if (present(m)) then
      basis%has_sym_m = .true.
      basis%m = m

#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, "<basis%has_sym_m> = ", basis%has_sym_m
      write (STDERR, *) DBG, "<basis%m> = ", basis%m
#endif
    else
      basis%has_sym_m = .false.
      basis%m = 0

#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, "<basis%has_sym_m> = ", basis%has_sym_m
      write (STDERR, *) DBG, "<basis%m> = ", basis%m
      write (STDERR, *) DBG, "arbitrarily set <basis%m>"
#endif
    end if

    ! if <parity> symmetry imposed, set <basis%has_sym_parity>, <basis%parity>
    ! else, set <basis%has_sym_parity>, and arbitrarily set <basis%parity>
    if (present(parity)) then
      basis%has_sym_parity = .true.
      basis%parity = parity

#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, &
          "<basis%has_sym_parity> = ", basis%has_sym_parity
      write (STDERR, *) DBG, "<basis%parity> = ", basis%parity
#endif
    else
      basis%has_sym_parity = .false.
      basis%parity = 0

#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, &
          "<basis%has_sym_parity> = ", basis%has_sym_parity
      write (STDERR, *) DBG, "<basis%parity> = ", basis%parity
      write (STDERR, *) DBG, "arbitrarily set <basis%parity>"
#endif
    end if

    ! if <l> symmetry imposed, set <basis%has_sym_l>, <basis%l>
    ! else, set <basis%has_sym_l>, and arbitrarily set <basis%l>
    if (present(l)) then
      basis%has_sym_l = .true.
      basis%l = l

#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, "<basis%has_sym_l> = ", basis%has_sym_l
      write (STDERR, *) DBG, "<basis%l> = ", basis%l
#endif
    else
      basis%has_sym_l = .false.
      basis%l = 0

#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, "<basis%has_sym_l> = ", basis%has_sym_l
      write (STDERR, *) DBG, "<basis%l> = ", basis%l
      write (STDERR, *) DBG, "arbitrarily set <basis%l>"
#endif
    end if

    ! set <l_max>
    basis%l_max = l_max

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "<basis%l_max> = ", basis%l_max
#endif

    ! allocate <n_basis_l>, <alpha_l>
    allocate(basis%n_basis_l(0:l_max))
    allocate(basis%alpha_l(0:l_max))

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "allocated <basis%n_basis_l>"
    write (STDERR, *) DBG, "allocated <basis%alpha_l>"
#endif

    ! calculate <n_basis>, ignoring <n_basis_l> where symmetries require
    basis%n_basis = 0

    do ll = 0, l_max

#if (DEBUG_BASIS >= 4)
      write (STDERR, *) DBG, "<l> = ", ll
#endif

      if ((basis%has_sym_m .and. (ll < abs(basis%m))) &
          .or. (basis%has_sym_parity .and. ((-1)**ll /= basis%parity))) then
        basis%n_basis_l(ll) = 0

#if (DEBUG_BASIS >= 4)
        write (STDERR, *) DBG, "ignoring due to basis symmetry"
#endif
      else
        basis%n_basis_l(ll) = n_basis_l(ll)

        if (basis%has_sym_m) then
          basis%n_basis = basis%n_basis + n_basis_l(ll)
        else
          basis%n_basis = basis%n_basis + (2*ll + 1)*n_basis_l(ll)
        end if
      end if
    end do

    ! set <alpha_l>
    basis%alpha_l(:) = alpha_l(:)

#if (DEBUG_BASIS >= 2)
    if (basis%has_sym_m) then
      write (STDERR, *) DBG, "<l>, <basis%n_basis_l>, <basis%alpha_l>"
      do ll = 0, basis%l_max
        write (STDERR, *) DBG, ll, basis%n_basis_l(ll), basis%alpha_l(ll)
      end do
    else
      write (STDERR, *) DBG, "<l>, <2*ll + 1>, <basis%n_basis_l>, <basis%alpha_l>"
      do ll = 0, basis%l_max
        write (STDERR, *) DBG, &
            ll, 2*ll + 1, basis%n_basis_l(ll), basis%alpha_l(ll)
      end do
    end if

    write (STDERR, *) DBG, "<basis%n_basis> = ", basis%n_basis
#endif

    ! if basis is empty, deallocate memory and return an error code
    if (basis%n_basis == 0) then
      deallocate(basis%n_basis_l)
      deallocate(basis%alpha_l)

      i_err = 1

#if (DEBUG_BASIS >= 1)
      write (STDERR, *) DBG, ERR, "basis is empty"
      write (STDERR, *) DBG, ERR, "<i_err> = ", i_err
      write (STDERR, *) DBG, ERR, "exiting subroutine setup_basis()"
#endif
      return
    end if

    ! allocate <k_list>, <l_list>, <m_list>
    allocate(basis%k_list(basis%n_basis))
    allocate(basis%l_list(basis%n_basis))
    allocate(basis%m_list(basis%n_basis))

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "allocated <basis%k_list>"
    write (STDERR, *) DBG, "allocated <basis%l_list>"
    write (STDERR, *) DBG, "allocated <basis%m_list>"
#endif

    ! set <k_list>, <l_list>, <m_list>
    ii = 0

    do mm = 0, l_max
#if (DEBUG_BASIS >= 4)
      write (STDERR, *) DBG, "<m> = ", mm
#endif

      if (basis%has_sym_m .and. (mm /= abs(basis%m))) then
#if (DEBUG_BASIS >= 4)
        write (STDERR, *) DBG, "ignoring due to symmetry"
        write (STDERR, *) DBG, "cycling"
#endif
        cycle
      end if

      do ll = 0, l_max
#if (DEBUG_BASIS >= 4)
        write (STDERR, *) DBG, "<l> = ", ll
#endif

        ! case of <m> = +<mm>
        if ((ll < abs(mm)) .or. (mm /= basis%m)) then
#if (DEBUG_BASIS >= 4)
          write (STDERR, *) DBG, "ignoring due to symmetry"
          write (STDERR, *) DBG, "cycling"
#endif
          cycle
        end if

        if (basis%n_basis_l(ll) >= 1) then
          do kk = 1, basis%n_basis_l(ll)
            basis%k_list(ii+kk) = kk
            basis%l_list(ii+kk) = ll
            basis%m_list(ii+kk) = mm
          end do
        end if
        ii = ii + basis%n_basis_l(ll)

        ! avoid duplication for case of <mm> = 0
        if (mm == 0) then
#if (DEBUG_BASIS >= 4)
          write (STDERR, *) DBG, "ignoring to avoid duplication"
          write (STDERR, *) DBG, "cycling"
#endif
          cycle
        end if

        ! case of <m> = -<mm>
        if ((ll < abs(mm)) .or. (mm /= -basis%m)) then
#if (DEBUG_BASIS >= 4)
          write (STDERR, *) DBG, "ignoring due to symmetry"
          write (STDERR, *) DBG, "cycling"
#endif
          cycle
        end if

        if (basis%n_basis_l(ll) >= 1) then
          do kk = 1, basis%n_basis_l(ll)
            basis%k_list(ii+kk) = kk
            basis%l_list(ii+kk) = ll
            basis%m_list(ii+kk) = -mm
          end do
        end if
        ii = ii + basis%n_basis_l(ll)
      end do
    end do

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, &
        "<i>, <basis%k_list>, <basis%l_list>, <basis%m_list>"
    do ii = 1, basis%n_basis
      write (STDERR, *) DBG, &
          ii, basis%k_list(ii), basis%l_list(ii), basis%m_list(ii)
    end do
#endif

    ! set <n_r> to zero, indicating radial basis functions have not yet been
    ! plotted on a grid
    basis%n_r = 0

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "<basis%n_r> = ", basis%n_r
#endif

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "end subroutine setup_states()"
#endif

  end subroutine setup_states

end module basis
