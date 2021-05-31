!>
module parameters

  ! debug parameters
  ! - <STDERR>: file unit for stderr output;
  ! - <DBG>: prefix every debug statement with this string;
  ! - <ERR>: prefix every error debug statement with this string;
  ! - <DEBUG>: verbosity of debug statements;
  !   - 0: none;
  !   - 1: control flow (entering and exiting subroutines/functions);
  !   - 2: allocation, scalar assignment, short array assignment, etc;
  !   - 3: local variable assignment, allocation, inspecting loops, etc;
  !   - 4: medium array assignment;
  !   - 5: longer array assignment.
#define STDERR 0
#define DBG "[debug] "
#define ERR "[error] "
#define DEBUG 2

  implicit none

  ! numeric parameters
  ! - <TOL>: double precision tolerance value.
  double precision , parameter :: TOL = 1.0D-10

end module parameters
