//
// notation:
// - S(n) = {0, ..., n-1}
// notes:
// - 0-indexing is used rather than 1-indexing

#include "group_symmetric.h"

#include <stdint.h>
#include <assert.h>

// integer type
typedef unsigned int u_int;


// permutation types

// transposition:
// - transpose `left -> right`, `right -> left`
// notes:
// - a transposition is simply a cycle with `length = 2`
// - the ordering of `left`, `right` is arbitrary, so we choose `left < right`
// example:
// - `Transposition t = {3, 9}` corresponds to the cycle `(3 9)`
struct Transposition
{
  u_int left;
  u_int right;
};
typedef struct Transposition Transposition;

// cycle:
// - cycle `length` unsigned integers, with `map[i] -> map[i+1]` for all
//   `i in S(length-1)` and `map[length] -> map[0]`
// constatins:
// - `map` must have `length` elements
// example:
// - `Cycle c = {3, {1, 5, 3}}` corresponds to the cycle `(1 5 3)` in cycle
//   notation
struct Cycle
{
  u_int length;
  u_int * map;
};
typedef struct Cycle Cycle;

// permutation:
// - permute `degree` unsigned integers, with `i -> map[i]` for all
//   `i in S(degree)`
// constraints:
// - `map` must have `length` elements
// - `map[i]` must be a bijective function, `map : S(degree) -> S(degree)`
// example:
// - `Permutation c = {4, {1, 4, 2, 3}}` corresponds to the permutation
//   `1 4 2 3`, in one-line notation, `(1 2 3 4 -> 1 4 2 3)` in two-line
//   notation, and `(2 4 3)` in cycle notation
struct Permutation
{
  u_int degree;
  u_int * map;
};
typedef struct Permutation Permutation;

// actions of permutations types

// transpose:
u_int transpose(Transposition * t, u_int i)
{
  if (i == t->left)
  {
    return (t->right);
  }

  if (i == t->right)
  {
    return (t->left);
  }

  return i;
}

// cycle:
u_int cycle(Cycle * c, u_int i)
{
  u_int idx = 0;

  // case of `i` being the not-last element of the cycle
  while (idx < c->length - 1)
  {
    if (i == c->map[idx])
    {
      return (c->map[idx+1]);
    }
    ++idx;
  }

  // case of `i` being the last element of the cycle
  if ((c->length > 0) && (i == c->map[c->length-1]))
  {
    return (c->map[0]);
  }

  // case of `i` not being in the cycle
  return i;
}

// permute:
u_int permute(Permutation * p, u_int i)
{
  assert(i < p->degree);

  return (p->map[i]);
}

// decompositions


// view


// constructors

// new_transposition:
// notes:
// - the ordering of `left`, `right` is arbitrary, so we choose `left < right`
new_transposition(Transposition * t, u_int left, u_int right)
{
  if (left <= right)
  {
    t->left = left;
    t->right = right;
  }
  else
  {
    t->left = right;
    t->right = left;
  }
}

// free
