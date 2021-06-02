//
// notation:
// - S(n) = {0, ..., n-1}
// notes:
// - 0-indexing is used rather than 1-indexing

#include "group_symmetric.h"

// permutation types

// transposition:
// - transpose `left -> right`, `right -> left`
// constatins:
// - `left`, `right` must be distinct
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

// cycle:
// - cycle `length` unsigned integers, with `map[i] -> map[i+1]` for all
//   `i in S(length-1)` and `map[length] -> map[0]`
// constatins:
// - `length` must be at least 2
// - `map` must have `length` elements
// example:
// - `Cycle c = {3, {1, 5, 3}}` corresponds to the cycle `(1 5 3)` in cycle
//   notation
struct Cycle
{
  u_int length;
  u_int * map;
};

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


// actions of permutations types

// transpose:
u_int transpose(Transposition const * const t, u_int const i)
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
u_int cycle(Cycle const * const c, u_int const i)
{
  u_int ii = 0;

  // case of `i` being the not-last element of the cycle
  while (ii < c->length - 1)
  {
    if (i == c->map[ii])
    {
      return (c->map[ii+1]);
    }
    ii++;
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
u_int permute(Permutation const * const p, u_int const i)
{
  assert(i < p->degree);

  return (p->map[i]);
}

// decompositions


// view

// view_transposition:
char * view_transposition(Transposition const * const t)
{
  assert(valid_transposition(t->left, t->right));

  // format string
  char const fmt[] = "(%i %i)";

  // safely determine size of `str` needed
  int const len = snprintf(NULL, 0, fmt, t->left, t->right);

  // construct string
  char * str = (char *)malloc(sizeof(char) * (size_t)(len + 1));
  int const mark = snprintf(str, (size_t)(len + 1), fmt, t->left, t->right);

  // ensure all characters were written
  assert((size_t)mark == strlen(str));

  return str;
}

// view_cycle:
char * view_cycle(Cycle const * const c)
{
  assert(valid_cycle(c->length, c->map));

  // format string
  char const fmt_l[] = "(";
  char const fmt_element[] = "%i";
  char const fmt_between[] = " ";
  char const fmt_r[] = ")";

  // safely determine size of `str` needed
  int len = snprintf(NULL, 0, fmt_l);
  for (u_int ii = 0; ii < c->length; ii++)
  {
    len += snprintf(NULL, 0, fmt_element, c->map[ii]);

    if (ii < c->length - 1)
    {
      len += snprintf(NULL, 0, fmt_between);
    }
  }
  len += snprintf(NULL, 0, fmt_r);

  // construct string
  char * str = (char *)malloc(sizeof(char) * (size_t)(len + 1));

  int mark = snprintf(str, (size_t)(len+1), fmt_l);
  for (u_int ii = 0; ii < c->length; ii++)
  {
    mark += snprintf(str+mark, (size_t)(len-mark+1), fmt_element, c->map[ii]);

    if (ii < c->length - 1)
    {
      mark += snprintf(str+mark, (size_t)(len-mark+1), fmt_between);
    }

  }
  mark += snprintf(str+mark, (size_t)(len-mark+1), fmt_r);

  // ensure all characters were written
  assert((size_t)mark == strlen(str));

  return str;
}

// view_permutation:
char * view_permutation(Permutation const * const p)
{
  assert(valid_permutation(p->degree, p->map));

  // format string
  char const fmt_l[] = "(";
  char const fmt_element[] = "(%i->%i)";
  char const fmt_between[] = " ";
  char const fmt_r[] = ")";

  // safely determine size of `str` needed
  int len = snprintf(NULL, 0, fmt_l);
  for (u_int ii = 0; ii < p->degree; ii++)
  {
    len += snprintf(NULL, 0, fmt_element, ii, p->map[ii]);

    if (ii < p->degree - 1)
    {
      len += snprintf(NULL, 0, fmt_between);
    }
  }
  len += snprintf(NULL, 0, fmt_r);

  // construct string
  char * str = (char *)malloc(sizeof(char) * (size_t)(len + 1));

  int mark = snprintf(str, (size_t)(len+1), fmt_l);
  for (u_int ii = 0; ii < p->degree; ii++)
  {
    mark += snprintf(str+mark, (size_t)(len-mark+1), fmt_element,
                     ii, p->map[ii]);

    if (ii < p->degree - 1)
    {
      mark += snprintf(str+mark, (size_t)(len-mark+1), fmt_between);
    }

  }
  mark += snprintf(str+mark, (size_t)(len-mark+1), fmt_r);

  // ensure all characters were written
  assert((size_t)mark == strlen(str));

  return str;
}

// constructors

// new_transposition:
// notes:
// - the ordering of `left`, `right` is arbitrary, so we choose `left < right`
void new_transposition(Transposition * const t, u_int const left,
                       u_int const right)
{
  // ensure `left`, `right` are distinct
  assert(valid_transposition(left, right));

  // order `t->left = min(left, right)`, `t->right = max(left, right)`
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

// new_cycle:
void new_cycle(Cycle * const c, u_int const length, u_int const * const map)
{
  // ensure `map` has at least 2 elements, and that all elements of `map` are
  // unique
  assert(valid_cycle(length, map));

  // set `c->length` to `length`
  c->length = length;

  // free `c->map` if it is allocated
  if (c->map != NULL) free(c->map);
  c->map = NULL;

  // copy`c->map` from `map` if `length` is non-zero
  if (length > 0)
  {
    size_t const n_bytes = sizeof(u_int) * (size_t)length;
    c->map = (u_int *)malloc(n_bytes);
    memcpy(c->map, map, n_bytes);
  }
}

// new_permutation:
void new_permutation(Permutation * const p, u_int const degree,
                     u_int const * const map)
{
  // ensure all elements of `map` are unique, and are in `S(degree)`
  assert(valid_permutation(degree, map));

  // set `p->degree` to `degree`
  p->degree = degree;

  // free `p->map` if it is allocated
  if (p->map != NULL) free(p->map);
  p->map = NULL;

  // copy `p->map` from `map` if `degree` is non-zero
  if (degree > 0)
  {
    size_t const n_bytes = sizeof(u_int) * (size_t)degree;
    p->map = (u_int *)malloc(n_bytes);
    memcpy(p->map, map, n_bytes);
  }
}

// free

// free_transposition:
// notes:
// - transpositions have no allocated memory to free
void free_transposition(Transposition * const t)
{
  t->left = 0;
  t->right = 0;
}

// free_cycle:
void free_cycle(Cycle * const c)
{
  // free `c->map` if it is allocated
  if (c->map != NULL) free(c->map);
  c->map = NULL;

  // set `c->length` to
  c->length = 0;
}

// free_permutation:
void free_permutation(Permutation * const p)
{
  // free `p->map` if it is allocated
  if (p->map != NULL) free(p->map);
  p->map = NULL;

  // set `p->degree`
  p->degree = 0;
}

// validators

// valid_transposition:
bool valid_transposition(u_int const left, u_int const right)
{
  // ensure `left`, `right` are distinct
  bool valid = (left != right);

  return valid;
}

// valid_cycle:
bool valid_cycle(u_int const length, u_int const * const map)
{
  bool valid = true;

  // ensure `map` has at least 2 elements
  valid = valid && (length >= 2);

  // ensure all elements of `map` are unique
  u_int ii = 0;

  while (valid && ii < length)
  {
    u_int jj = ii + 1;

    while (valid && jj < length)
    {
      // ensure `map[ii]` is unique in `map`
      valid = valid && (map[ii] != map[jj]);
      jj++;
    }

    ii++;
  }

  return valid;
}

// valid_permutation:
bool valid_permutation(u_int const degree, u_int const * const map)
{
  // ensure all elements of `map` are unique, and are in `S(degree)`
  bool valid = true;
  u_int ii = 0;

  while (valid && ii < degree)
  {
    u_int jj = ii + 1;

    // ensure `map[ii] in S(degree)`
    valid = valid && (map[ii] < degree);

    while (valid && jj < degree)
    {
      // ensure `map[ii]` is unique in `map`
      valid = valid && (map[ii] != map[jj]);
      jj++;
    }

    ii++;
  }

  return valid;
}
