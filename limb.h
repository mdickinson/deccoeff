#include <stdint.h>
#include <stdbool.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "Python.h"
#include "longintrepr.h"

typedef int32_t limb_t;
#define LIMB_ZERO ((limb_t)0)
#define LIMB_ONE ((limb_t)1)
#define LIMB_MAX ((limb_t)999999999)

#define LIMB_DIGITS (9)

/* definitions used for conversion from binary to decimal and back */

/* digitpair is an integer type capable of holding a pair of PyLong
   digits (i.e., it should hold any integer in the range [0, 2**30).
   digitpair_limb_t holds any integer in the range [0,
   2**30*LIMB_BASE); these types are used for decimal<->binary base
   conversions. */

typedef int32_t digitpair;
#define DIGIT_PAIR(a, b) (((digitpair)(a) << PyLong_SHIFT) | (b))
#define DIGIT_PAIR_SHIFT (2*PyLong_SHIFT)
#define DIGIT_PAIR_BASE ((digitpair)1 << DIGIT_PAIR_SHIFT)
#define DIGIT_PAIR_MASK (DIGIT_PAIR_BASE - 1)

/* arithmetic operations on limbs */

/* add with carry */
bool limb_adc(limb_t *, limb_t, limb_t, bool);
/* subtract with borrow */
bool limb_sbb(limb_t *, limb_t, limb_t, bool);
/* multiply and add, with two addends */
limb_t limb_fmaa(limb_t *, limb_t, limb_t, limb_t, limb_t);
/* divide, returning quotient and remainder */
limb_t limb_div(limb_t *, limb_t, limb_t, limb_t);

/* comparisons */
bool limb_eq(limb_t, limb_t);
bool limb_le(limb_t, limb_t);
bool limb_lt(limb_t, limb_t);

/* shifts, rotates, etc. */
/* shift right */
limb_t limb_sar(limb_t, Py_ssize_t);
/* shift left */
limb_t limb_sal(limb_t, Py_ssize_t);
/* rotate right and split */
limb_t limb_split(limb_t *, limb_t, Py_ssize_t);
/* rotate left and split */
limb_t limb_splitl(limb_t *, limb_t, Py_ssize_t);
/* select low digits */
limb_t limb_low(limb_t, Py_ssize_t);
/* index of most significant nonzero digit */
Py_ssize_t limb_dsr(limb_t);

/* helpers for base conversion

   given positive integers M and N, any x in range(0, M*N) can be
   represented uniquely in the form

     a*M + b,  0 <= a < N,  0 <= b < M

   and also in the form

     c*N + d, 0 <= c < M,  0 <= d < N.

  The following two functions convert from one representation to the
  other, in the particular case where M is 2**30 and N is LIMB_BASE.
  These two operations are exactly the primitive operations needed for
  binary<->decimal base conversion. */

digitpair limb_digitpair_swap(limb_t *, limb_t, digitpair);
limb_t digitpair_limb_swap(digitpair *, digitpair, limb_t);

/* functions for conversion to and from strings */
/* get a particular digit, as a character */
char limb_getdigit(limb_t, Py_ssize_t);
/* set a particular digit */
limb_t limb_setdigit(limb_t, Py_ssize_t, char);

/* hash of a single limb;  used for making deccoeff hashable */
unsigned long limb_hash(limb_t);

Py_ssize_t limbsize_from_longsize(Py_ssize_t);
Py_ssize_t longsize_from_limbsize(Py_ssize_t);
