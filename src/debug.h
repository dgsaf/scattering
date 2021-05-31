#ifndef DEBUG_H
#define DEBUG_H

/*
  debug parameters

  <STDERR>: file unit for stderr output;
  <DBG>: prefix every debug statement with this string;
  <ERR>: prefix every error debug statement with this string;
  <DEBUG>: verbosity of debug statements;
  - 0: none;
  - 1: control flow (entering and exiting subroutines/functions);
  - 2: allocation, scalar assignment, short array assignment, etc;
  - 3: local variable assignment, allocation, inspecting loops, etc;
  - 4: medium array assignment;
  - 5: longer array assignment.
*/

#define STDERR 0
#define DBG "[debug] "
#define ERR "[error] "
#define DEBUG 2

#endif
