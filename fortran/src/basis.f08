!>
module basis

#include "debug.h"

  ! debug compilation
  ! - <DEBUG_BASIS>: verbosity of debug statements, local to this file.
#define DEBUG_BASIS 3

  ! local debug verbosity cannot be lower than the global debug verbosity
#if ((defined DEBUG) && (defined DEBUG_BASIS) && (DEBUG > DEBUG_BASIS))
#undef DEBUG_BASIS
#define DEBUG_BASIS DEBUG
#endif

  use parameters

  implicit none

  private
  public setup_states, setup_radial, valid_states, &
      partial_waves, overlap, kinetic

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
    double precision , allocatable :: radial(:, :, :)
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

#if (DEBUG_BASIS >= 3)
      write (STDERR, *) DBG, "<l> = ", ll
#endif

      if ((basis%has_sym_m .and. (ll < abs(basis%m))) &
          .or. (basis%has_sym_parity .and. ((-1)**ll /= basis%parity))) then
        basis%n_basis_l(ll) = 0

#if (DEBUG_BASIS >= 3)
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
#if (DEBUG_BASIS >= 3)
      write (STDERR, *) DBG, "<m> = ", mm
#endif

      if (basis%has_sym_m .and. (mm /= abs(basis%m))) then
#if (DEBUG_BASIS >= 3)
        write (STDERR, *) DBG, "ignoring due to symmetry"
        write (STDERR, *) DBG, "cycling"
#endif
        cycle
      end if

      do ll = 0, l_max
#if (DEBUG_BASIS >= 3)
        write (STDERR, *) DBG, "<l> = ", ll
#endif

        ! case of <m> = +<mm>
        if ((ll < abs(mm)) .or. (mm /= basis%m)) then
#if (DEBUG_BASIS >= 3)
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
#if (DEBUG_BASIS >= 3)
          write (STDERR, *) DBG, "ignoring to avoid duplication"
          write (STDERR, *) DBG, "cycling"
#endif
          cycle
        end if

        ! case of <m> = -<mm>
        if ((ll < abs(mm)) .or. (mm /= -basis%m)) then
#if (DEBUG_BASIS >= 3)
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

  ! setup_radial
  !
  ! For given <basis>, <n_r>, <r_grid>, calculates the radial functions
  ! > varphi_{k_{i}, l_{i}}(r) for i = 1, .., <n_basis>
  ! on the radial values specified in the grid.
  !
  ! Requires the following variables in <basis> to have already been setup:
  ! - <l_max>;
  ! - <n_basis_l>;
  ! - <alpha_l>;
  ! - <n_basis>.
  ! That is, a call to setup_states() should have already been made.
  !
  ! Returns an error code, <i_err>, where:
  ! - 0 indicates successful execution;
  ! - 1 indicates invalid arguments.
  subroutine setup_radial (basis, n_r, r_grid, i_err)
    type(t_basis) , intent(inout) :: basis
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    integer , intent(out) :: i_err
    double precision , allocatable :: norm(:)
    double precision :: alpha_grid(n_r)
    double precision :: alpha
    integer :: n_b_l
    integer :: ii, kk, ll

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "subroutine setup_radial()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (.not. valid_states(basis)) then
      i_err = 1
    end if

    if (n_r < 1) then
      i_err = 1
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<n_r> < 1"
#endif
    end if

    ! handle invalid arguments
    if (i_err /= 0) then
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_BASIS >= 1)
      write (STDERR, *) DBG, ERR, "arguments are invalid"
      write (STDERR, *) DBG, ERR, &
          "<basis> array variables will be left un-allocated"
      write (STDERR, *) DBG, ERR, "exiting subroutine setup_radial()"
#endif
      return
    end if

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "arguments are valid"
#endif

    ! set <n_r>
    basis%n_r = n_r

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "<basis%n_r> = ", basis%n_r
#endif

    ! allocate <r_grid>, <radial>
    allocate(basis%r_grid(n_r))
    allocate(basis%radial(n_r, maxval(basis%n_basis_l(:)), 0:basis%l_max))

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "allocated <basis%r_grid>"
    write (STDERR, *) DBG, "allocated <basis%radial>"
#endif

    ! basis: set <r_grid>
    basis%r_grid(:) = r_grid(:)

#if (DEBUG_BASIS >= 5)
    write (STDERR, *) DBG, "<radial index>, <basis%r_grid>"
    do ii = 1, basis%n_r
      write (STDERR, *) DBG, ii, basis%r_grid(ii)
    end do
#endif

    ! loop over <l>, basis: set <radial>
    basis%radial(:, :, :) = 0.0d0

    do ll = 0, basis%l_max

#if (DEBUG_BASIS >= 3)
      write (STDERR, *) DBG, "<l> = ", ll
#endif

      ! in-line <n_basis_l>, <alpha> for current <l>
      n_b_l = basis%n_basis_l(ll)
      alpha = basis%alpha_l(ll)

#if (DEBUG_BASIS >= 3)
      write (STDERR, *) DBG, "<n_b_l> = ", n_b_l
      write (STDERR, *) DBG, "<alpha> = ", alpha
#endif

      ! if there are no basis states for this value of <l>, then cycle
      if (n_b_l == 0) then

#if (DEBUG_BASIS >= 3)
        write (STDERR, *) DBG, "no basis states for this value of <l>"
        write (STDERR, *) DBG, "cycling"
#endif
        cycle
      end if

      ! allocate normalisation constant array
      allocate(norm(n_b_l))

#if (DEBUG_BASIS >= 3)
      write (STDERR, *) DBG, "allocated <norm>"
#endif

      ! recurrence relation for basis normalisation constants
#if (DEBUG_BASIS >= 3)
      write (STDERR, *) DBG, "<k>, <norm>"
#endif

      if (n_b_l >= 1) then
        norm(1) = sqrt(alpha / dble((ll+1) * gamma(dble((2*ll)+2))))

#if (DEBUG_BASIS >= 3)
        write (STDERR, *) DBG, 1, norm(1)
#endif
      end if

      if (n_b_l >= 2) then
        do kk = 2, n_b_l
          norm(kk) = norm(kk-1) * sqrt(dble((kk-1) * (kk-1+ll)) &
              / dble((kk+ll) * (kk+(2*ll))))

#if (DEBUG_BASIS >= 3)
          write (STDERR, *) DBG, kk, norm(kk)
#endif
        end do
      end if

      ! in-line <alpha_grid> for current <l>
      alpha_grid(:) = alpha * r_grid(:)

#if (DEBUG_BASIS >= 5)
      write (STDERR, *) DBG, "<radial index>, <alpha_grid>"
      do ii = 1, basis%n_r
        write (STDERR, *) DBG, ii, alpha_grid(ii)
      end do
#endif

      ! recurrence relation for basis functions
      if (n_b_l >= 1) then
        basis%radial(:, 1, ll) = ((2.0d0 * alpha_grid(:)) ** (ll+1)) &
            * exp(-alpha_grid(:))
      end if

      if (n_b_l >= 2) then
        basis%radial(:, 2, ll) = 2.0d0 * (dble(ll+1) - alpha_grid(:)) &
            * basis%radial(:, 1, ll)
      end if

      if (n_b_l >= 3) then
        do kk = 3, n_b_l
          basis%radial(:, kk, ll) = &
              ((2.0d0 * (dble(kk-1+ll) - alpha_grid(:)) &
              * basis%radial(:, kk-1, ll)) &
              - dble(kk+(2*ll)-1) * basis%radial(:, kk-2, ll)) &
              / dble(kk-1)
        end do
      end if

#if (DEBUG_BASIS >= 5)
      write (STDERR, *) DBG, "<radial index>, <basis%radial> (un-normalised)"
      do ii = 1, basis%n_r
        write (STDERR, *) DBG, ii, basis%radial(ii, :, ll)
      end do
#endif

      ! scale basis radial functions by normalisation constants
      if (n_b_l >= 1) then
        do kk = 1, n_b_l
          basis%radial(:, kk, ll) = basis%radial(:, kk, ll) * norm(kk)
        end do
      end if

#if (DEBUG_BASIS >= 5)
      write (STDERR, *) DBG, "<radial index>, <basis%radial>"
      do ii = 1, basis%n_r
        write (STDERR, *) DBG, ii, basis%radial(ii, :, ll)
      end do
#endif

      ! deallocate norm
      deallocate(norm)

#if (DEBUG_BASIS >= 3)
      write (STDERR, *) DBG, "deallocated <norm>"
#endif

    end do

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "end subroutine setup_radial()"
#endif

  end subroutine setup_radial

  ! valid_states
  !
  ! For given <basis>, determines if the basis variables are valid. A <basis> is
  ! invalid if:
  ! - abs(<basis%parity>) /= 1;
  ! - <basis%l_max> < abs(<basis%m>);
  ! - <basis%l_max> < 0;
  ! - <basis%n_basis> < 1;
  ! - <basis%l_max> < 0;
  ! - any(<basis%n_basis_l(:)> < 0);
  ! - any(<basis%alpha_l(:)> < TOL).
  logical function valid_states (basis)
    type(t_basis) , intent(in) :: basis

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "function valid_states()"
#endif

    ! check if <basis> is valid
    valid_states = .true.

    if (basis%n_basis < 1) then
      valid_states = .false.
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<basis%n_basis> < 1"
#endif
    end if

    if (basis%has_sym_parity) then
      if (abs(basis%parity) /= 1) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "abs(<basis%parity>) /= 1"
#endif
      end if

      if (any(((-1) ** basis%l_list(:)) /= basis%parity)) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, &
            "any(((-1) ** <basis%l_list(:)>) /= <basis%parity>)"
#endif
      end if
    end if

    if (basis%has_sym_m) then
      if (basis%l_max < abs(basis%m)) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "<basis%l_max> < abs(<basis%m>)"
#endif
      end if

      if (any(basis%m_list(:) /= basis%m)) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "any(<basis%m_list(:)> /= <basis%m>)"
#endif
      end if
    end if

    if (basis%has_sym_l) then
      if (basis%l < 0) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "<basis%l> < 0"
#endif
      end if

      if (basis%l > basis%l_max) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "<basis%l> > <basis%l_max>"
#endif
      end if

      if (any(basis%l_list(:) /= basis%l)) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "any(<basis%l_list(:)> /= <basis%l>)"
#endif
      end if

      if (basis%has_sym_m .and. (basis%l < abs(basis%m))) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "<basis%l> < abs(<basis%m>)"
#endif
      end if

      if (basis%has_sym_parity .and. (((-1) ** basis%l) /= basis%parity)) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "(-1) ** <basis%l> /= <basis%parity>"
#endif
      end if
    end if

    if (basis%l_max < 0) then
      valid_states = .false.
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<basis%l_max> < 0"
#endif
    else
      if (any(basis%n_basis_l(:) < 0)) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "any(<basis%n_basis_l(:)> < 0)"
#endif
      end if

      if (any(basis%alpha_l(:) < TOL)) then
        valid_states = .false.
#if (DEBUG_BASIS >= 2)
        write (STDERR, *) DBG, ERR, "any(<basis%alpha_l(:)> < TOL)"
#endif
      end if
    end if

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "end function valid_states()"
#endif

    return
  end function valid_states

  ! partial_waves
  !
  ! For given <basis>, <C>, where <C> is a coefficient matrix for a set of
  ! states, in terms of <basis>, and where the <basis> has an imposed <m>
  ! symmetry, calculate the partial waves of the states.
  !
  ! Returns an error code, <i_err>, where:
  ! - 0 indicates successful execution;
  ! - 1 indicates invalid arguments.
  subroutine partial_waves (basis, C, pw, i_err)
    type(t_basis) , intent(in) :: basis
    double precision , intent(in) :: C(basis%n_basis, basis%n_basis)
    double precision , intent(out) :: pw(basis%n_r, basis%l_max, basis%n_basis)
    integer , intent(out) :: i_err
    integer :: ii, jj, nn, kk, ll

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "subroutine partial_waves()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (.not. valid_states(basis)) then
      i_err = 1
    end if

    if (basis%n_r < 1) then
      i_err = 1
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<basis%n_r> < 1"
#endif
    end if

    if (.not. basis%has_sym_m) then
      i_err = 1
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<basis%has_sym_m> == .false."
#endif
    end if

    ! handle invalid arguments
    if (i_err /= 0) then
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_BASIS >= 1)
      write (STDERR, *) DBG, ERR, "arguments are invalid"
      write (STDERR, *) DBG, ERR, "exiting subroutine partial_waves()"
#endif
      return
    end if

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "arguments are valid"
#endif

    ! initialise <pw>
    pw(:, :, :) = 0.0d0

    ! calculate radial functions of expanded states
    do ii = 1, basis%n_basis

#if (DEBUG_BASIS >= 5)
      write (STDERR, *) DBG, "<i> = " ii
#endif

#if (DEBUG_BASIS >= 5)
      write (STDERR, *) DBG, "<radial index>, <pw>"
#endif

      do jj = 1, basis%n_r
        do nn = 1, basis%n_basis
          kk = basis%k_list(nn)
          ll = basis%l_list(nn)

          pw(jj, ll, ii) = pw(jj, ll, ii) &
              + (C(nn, ii) * basis%radial(jj, kk, ll))
        end do

#if (DEBUG_BASIS >= 5)
        write (STDERR, *) DBG, jj, pw(jj, :, ii)
#endif

      end do
    end do

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "end subroutine partial_waves()"
#endif

  end subroutine partial_waves

  ! overlap
  !
  ! For given <basis>, calculate the overlap matrix elements
  ! > B_{i, j} = < phi_{i} | phi_{j} > for i, j = 1, .., <n_basis>
  ! using analytic properties of the Laguerre basis.
  !
  ! Returns an error code, <i_err>, where:
  ! - 0 indicates successful execution;
  ! - 1 indicates invalid arguments.
  subroutine overlap (basis, B, i_err)
    type(t_basis) , intent(in) :: basis
    double precision , intent(out) :: B(basis%n_basis, basis%n_basis)
    integer , intent(out) :: i_err
    integer :: k_i, l_i, m_i, k_j, l_j, m_j
    integer :: ii, jj

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "subroutine overlap()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (.not. valid_states(basis)) then
      i_err = 1
    end if

    ! handle invalid arguments
    if (i_err /= 0) then
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_BASIS >= 1)
      write (STDERR, *) DBG, ERR, "arguments are invalid"
      write (STDERR, *) DBG, ERR, "exiting subroutine overlap()"
#endif
      return
    end if

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "arguments are valid"
#endif

    ! initialise <B>
    B(:, :) = 0.0d0

    ! calculate overlap matrix elements
#if (DEBUG_BASIS >= 4)
    write (STDERR, *) DBG, "<i>, <j>, <B(i, j)>"
#endif

    do jj = 1, basis%n_basis
      k_j = basis%k_list(jj)
      l_j = basis%l_list(jj)
      m_j = basis%m_list(jj)

      do ii = 1, basis%n_basis
        k_i = basis%k_list(ii)
        l_i = basis%l_list(ii)
        m_i = basis%m_list(ii)

        if ((l_i == l_j) .and. (m_i == m_j)) then
          if (k_i == k_j) then
            B(ii, jj) = 1.0d0
          else if (k_i + 1 == k_j) then
            B(ii, jj) = - 0.5d0 * sqrt(1.0d0 - &
                (dble(l_i * (l_i + 1)) / dble((k_i + l_i) * (k_i + l_i + 1))))
          else if (k_i == k_j + 1) then
            B(ii, jj) = - 0.5d0 * sqrt(1.0d0 - &
                (dble(l_j * (l_j + 1)) / dble((k_j + l_j) * (k_j + l_j + 1))))
          end if

#if (DEBUG_BASIS >= 4)
          write (STDERR, *) DBG, ii, jj, B(ii, jj)
#endif
        end if
      end do
    end do

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "end subroutine overlap()"
#endif

  end subroutine overlap

  ! kinetic
  !
  ! For given <basis>, calculate the kinetic matrix elements
  ! > K_{i, j} = < phi_{i} | K | phi_{j} > for i, j = 1, .., <n_basis>
  ! using analytic properties of the Laguerre basis.
  !
  ! Returns an error code, <i_err>, where:
  ! - 0 indicates successful execution;
  ! - 1 indicates invalid arguments.
  subroutine kinetic (basis, K, i_err)
    type(t_basis) , intent(in) :: basis
    double precision , intent(out) :: K(basis%n_basis, basis%n_basis)
    integer , intent(out) :: i_err
    double precision :: alpha
    integer :: k_i, l_i, m_i, k_j, l_j, m_j
    integer :: ii, jj

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "subroutine kinetic()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (.not. valid_states(basis)) then
      i_err = 1
    end if

    ! handle invalid arguments
    if (i_err /= 0) then
#if (DEBUG_BASIS >= 2)
      write (STDERR, *) DBG, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_BASIS >= 1)
      write (STDERR, *) DBG, ERR, "arguments are invalid"
      write (STDERR, *) DBG, ERR, "exiting subroutine kinetic()"
#endif
      return
    end if

#if (DEBUG_BASIS >= 2)
    write (STDERR, *) DBG, "arguments are valid"
#endif

    ! initialise <K>
    K(:, :) = 0.0d0

    ! calculate kinetic matrix elements
#if (DEBUG_BASIS >= 4)
    write (STDERR, *) DBG, "<i>, <j>, <K(i, j)>"
#endif

    do jj = 1, basis%n_basis
      k_j = basis%k_list(jj)
      l_j = basis%l_list(jj)
      m_j = basis%m_list(jj)

      alpha = basis%alpha_l(l_j)

      do ii = 1, basis%n_basis
        k_i = basis%k_list(ii)
        l_i = basis%l_list(ii)
        m_i = basis%m_list(ii)

        if ((l_i == l_j) .and. (m_i == m_j)) then
          if (k_i == k_j) then
            K(ii, jj) = 0.5d0 * (alpha ** 2)
          else if (k_i + 1 == k_j) then
            K(ii, jj) = (alpha ** 2) * 0.25d0 * sqrt(1.0d0 - &
                (dble(l_i * (l_i + 1)) / dble((k_i + l_i) * (k_i + l_i + 1))))
          else if (k_i == k_j + 1) then
            K(ii, jj) = (alpha ** 2) * 0.25d0 * sqrt(1.0d0 - &
                (dble(l_j * (l_j + 1)) / dble((k_j + l_j) * (k_j + l_j + 1))))
          end if

#if (DEBUG_BASIS >= 4)
          write (STDERR, *) DBG, ii, jj, K(ii, jj)
#endif
        end if
      end do
    end do

#if (DEBUG_BASIS >= 1)
    write (STDERR, *) DBG, "end subroutine kinetic()"
#endif

  end subroutine kinetic

end module basis
