//

#include "group_symmetric.h"

int main(int argc, char **argv)
{
  u_int n = 5;
  u_int c_map[] = {2, 5, 1, 3, 4};
  u_int p_map[] = {3, 5, 4, 2, 1};

  Transposition * t = (Transposition *)malloc(sizeof(t));
  new_transposition(t, 1, 3);
  printf("transposition: %s\n", view_transposition(t));
  for (u_int ii = 1; ii <= n; ii++)
  {
    printf("  %i -> %i\n", ii, transpose(t, ii));
  }

  Cycle * c = (Cycle *)malloc(sizeof(c));
  new_cycle(c, n, c_map);
  printf("cycle: %s\n", view_cycle(c));
  for (u_int ii = 1; ii <= n; ii++)
  {
    printf("  %i -> %i\n", ii, cycle(c, ii));
  }

  Permutation * p = (Permutation *)malloc(sizeof(p));
  new_permutation(p, n, p_map);
  printf("permutation: %s\n", view_permutation(p));
  for (u_int ii = 1; ii <= n; ii++)
  {
    printf("  %i -> %i\n", ii, permute(p, ii));
  }

  return 0;
}
