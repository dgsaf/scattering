//
// notation:
// - S(n) = {1, ..., n}
// notes:
// - 1-indexing is used for the outwards interface, rather than 0-indexing

#ifndef GROUP_SYMMETRIC_H
#define GROUP_SYMMETRIC_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

// unsigned integer type
typedef unsigned int u_int;

// permutation types
typedef struct Transposition Transposition;
typedef struct Cycle Cycle;
typedef struct Permutation Permutation;

// function declarations

// action
u_int transpose(Transposition const * const t, u_int const i);
u_int cycle(Cycle const * const c, u_int const i);
u_int permute(Permutation const * const p, u_int const i);

// view
char * view_transposition(Transposition const * const t);
char * view_cycle(Cycle const * const c);
char * view_permutation(Permutation const * const p);

// new
void new_transposition(Transposition * const t, u_int const left,
                       u_int const right);
void new_cycle(Cycle * const c, u_int const length, u_int const * const map);
void new_permutation(Permutation * const p, u_int const degree,
                     u_int const * const map);

// free
void free_transposition(Transposition * const t);
void free_cycle(Cycle * const c);
void free_permutation(Permutation * const p);

// validators
bool valid_transposition(u_int const left, u_int const right);
bool valid_cycle(u_int const length, u_int const * const map);
bool valid_permutation(u_int const degree, u_int const * const map);

#endif
