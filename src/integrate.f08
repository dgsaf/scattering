!>
module integrate

#include "debug.h"

  ! debug compilation
  ! - <DEBUG_INTEGRATE>: verbosity of debug statements, local to this file.
#define DEBUG_INTEGRATE 1

  ! local debug verbosity cannot be lower than the global debug verbosity
#if ((defined DEBUG) && (defined DEBUG_INTEGRATE) && (DEBUG > DEBUG_INTEGRATE))
#undef DEBUG_INTEGRATE
#define DEBUG_INTEGRATE DEBUG
#endif

  use parameters

  implicit none

contains

  function integrate_trapezoid (n, x, f) result (integral)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n), f(n)
    double precision :: w(n)
    double precision :: integral
    integer :: ii
    integer :: i_err

#if (DEBUG_INTEGRATE >= 1)
    write (STDERR, *) DBG, "function integrate_trapezoid()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (n < 2) then
      i_err = 1
#if (DEBUG_INTEGRATE >= 2)
      write (STDERR, *) DBG, ERR, "<n> < 2"
#endif
    end if

    ! handle invalid arguments
    if (i_err /= 0) then
#if (DEBUG_INTEGRATE >= 2)
      write (STDERR, *) DBG, ERR, "<i_err> = ", i_err
#endif

      integral = 0.0d0

#if (DEBUG_INTEGRATE >= 1)
      write (STDERR, *) DBG, ERR, "arguments are invalid"
      write (STDERR, *) DBG, ERR, "returning <integral> = 0"
      write (STDERR, *) DBG, ERR, "exiting function integrate_trapezoid()"
#endif
      return
    end if

#if (DEBUG_INTEGRATE >= 2)
    write (STDERR, *) DBG, "arguments are valid"
#endif

#if (DEBUG_INTEGRATE >= 5)
    write (STDERR, *) DBG, "<index>, <w>"
#endif

    ! calculate weights
    w(1) = (x(2) - x(1)) / 2.0d0

#if (DEBUG_INTEGRATE >= 5)
    write (STDERR, *) DBG, 1, w(1)
#endif

    do ii = 2, n-1
      w(ii) = (x(ii+1) - x(ii-1)) / 2.0d0

#if (DEBUG_INTEGRATE >= 5)
      write (STDERR, *) DBG, ii, w(ii)
#endif

    end do

    w(n) = (x(n) - x(n-1)) / 2.0d0

#if (DEBUG_INTEGRATE >= 5)
    write (STDERR, *) DBG, n, w(n)
#endif

    ! calculate integral
    integral = 0.0d0

    do ii = 1, n
      integral = integral + (w(ii) * f(ii))
    end do

#if (DEBUG_INTEGRATE >= 2)
    write (STDERR, *) DBG, "<integral> = ", integral
#endif

#if (DEBUG_INTEGRATE >= 1)
    write (STDERR, *) DBG, "end function integrate_trapezoid()"
#endif
    return

  end function integrate_trapezoid

end module integrate
