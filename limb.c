/*
   Base operations for Deccoeff type.  Definition of single limbs, and
   operations to work with single limbs.

   The two files limbs.c and limbs.h are supposed to encapsulate everything
   related to the representation of limbs; other code should be independent of
   the limb representation, only needing to know the number of digits in each
   limb.

   So for example it should be easy to change the implementation to use a
   16-bit limb with 4 digits per limb, or a 64-bit limb containing 18 or 19
   digits per limb (this might work well on 64-bit architectures like x86-64),
   or a binary-coded decimal or densely-packed decimal representation, or even
   a string of digits, or ...

   In particular, other code should not assume that limbs can be operated on
   using the usual arithmetic operators, or that they are comparable with < or
   == (some encodings may not be monotonic, or may have redundant encodings of
   the same integer, or may not even be encoded as a C integer type).
*/

#include <assert.h>
#include "limb.h"

/* for splint checking, we define Py_ssize_t, PY_SSIZE_T_MAX
   using explicit defines. */

#define LIMB_BASE (LIMB_MAX+1)

typedef int64_t double_limb_t;

static limb_t powers_of_ten[LIMB_DIGITS+1] = {
	1,
	10,
	100,
	1000,
	10000,
	100000,
	1000000,
	10000000,
	100000000,
	1000000000
};

/* XXX replace limb_error with Py_FataError eventually */

void
limb_error(const char *msg)
{
	fprintf(stderr, "%s\n", msg);
	abort();
}

/* add with carry: compute a + b + c, put result in *r and return
   the new carry. */

extern bool
limb_adc(limb_t *r, limb_t a, limb_t b, bool c)
{
	limb_t sum, test;
	sum = a + b + (c ? LIMB_ONE : LIMB_ZERO);
	test = sum - LIMB_BASE;
	if (test >= 0) {
		*r = test;
		return true;
	}
	else {
		*r = sum;
		return false;
	}
}

/* subtract with borrow:  compute a - (b + c); put result in *r
   and return the carry. */

extern bool
limb_sbb(limb_t *r, limb_t a, limb_t b, bool c)
{
	limb_t diff;
	diff = a - b - (c ? LIMB_ONE : LIMB_ZERO);
	if (diff < 0) {
		*r = diff + LIMB_BASE;
		return true;
	}
	else {
		*r = diff;
		return false;
	}
}

/* comparisons between limbs */

extern bool
limb_eq(limb_t a, limb_t b)
{
	return a == b;
}

extern bool
limb_le(limb_t a, limb_t b)
{
	return a <= b;
}

extern bool
limb_lt(limb_t a, limb_t b)
{
	return a < b;
}

/* multiply and add with two addends; this is useful for long
   multiplication.  Store the low part of the result in *low, and return
   the high part. */

extern limb_t
limb_fmaa(limb_t *low, limb_t a, limb_t b, limb_t c, limb_t d) {
	double_limb_t hilo;
	hilo = (double_limb_t)a * b + c + d;
	*low = (limb_t)(hilo%LIMB_BASE);
	return (limb_t)(hilo/LIMB_BASE);
}

/* division: divide high*BASE+low by c, giving a quotient (returned)
   and a remainder (stored in *rem).  Requires that high < c (and
   hence that c is nonzero). */

extern limb_t
limb_div(limb_t *rem, limb_t high, limb_t low, limb_t c) {
	double_limb_t hilo;
	if (high >= c)
		limb_error("invalid division");
	hilo = (double_limb_t)high * LIMB_BASE + low;
	*rem = (limb_t)(hilo%c);
	return (limb_t)(hilo/c);
}

/* return index of most significant digit of given limb;
   result is undefined if the limb is 0 */

extern Py_ssize_t
limb_dsr(limb_t x) {
	Py_ssize_t i;
	if (x == 0)
		limb_error("0 has no most significant digit");
	for (i=0; powers_of_ten[i] <= x; i++);
	return i;
}

/* shift a limb right n places (0 <= n <= LIMB_DIGITS) */

extern limb_t
limb_sar(limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_sar");
	return x / powers_of_ten[n];
}

/* rotate right and split; like limb_sar, but also puts piece that was shifted
   out into the top of *res */

extern limb_t
limb_split(limb_t *res, limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_split");
	*res = x % powers_of_ten[n] * powers_of_ten[LIMB_DIGITS-n];
	return x / powers_of_ten[n];
}

/* rotate left and split */

extern limb_t
limb_splitl(limb_t *res, limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_splitl");
	*res = x / powers_of_ten[LIMB_DIGITS-n];
	return x % powers_of_ten[LIMB_DIGITS-n] * powers_of_ten[n];
}

/* shift a limb left n places (0 <= n <= LIMB_DIGITS) */

extern limb_t
limb_sal(limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_sal");
	return x % powers_of_ten[LIMB_DIGITS-n] * powers_of_ten[n];
}

/* select the bottom n digits of a limb (0 <= n <= LIMB_DIGITS) */

extern limb_t
limb_low(limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_low");
	return x % powers_of_ten[n];
}

extern char
limb_getdigit(limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n < LIMB_DIGITS))
		limb_error("invalid digit position in limb_getdigit");
	return '0' + (char)((x / powers_of_ten[n]) % 10);
}

extern limb_t
limb_setdigit(limb_t x, Py_ssize_t n, char d) {
	if (!(0 <= n && n < LIMB_DIGITS))
		limb_error("invalid digit position in limb_setdigit");
	if (!('0' <= d && d <= '9'))
		limb_error("invalid digit in limb_setdigit");
	return x + powers_of_ten[n] * (limb_t)(d - limb_getdigit(x, n));
}

typedef int64_t digitpair_limb_t;

/* given limb a and digit pair b, write a*2**30 + b in the form
   c*LIMB_BASE + d, with c a digit pair and d a limb.  Return c and
   put d in *low. */

extern digitpair
limb_digitpair_swap(limb_t *low, limb_t a, digitpair b)
{
	digitpair_limb_t hilo;
	hilo = ((digitpair_limb_t)a << DIGIT_PAIR_SHIFT) + b;
	*low = (limb_t)(hilo%LIMB_BASE);
	return (digitpair)(hilo/LIMB_BASE);
}

/* reverse of the above: given a digit pair a and limb b, write LIMB_BASE*a +
   b in the form c*(2**30) + d, c a limb and d a digit pair.  Return c and put
   d in *low. */

extern limb_t
digitpair_limb_swap(digitpair *low, digitpair a, limb_t b)
{
	digitpair_limb_t hilo;
	hilo = (digitpair_limb_t)a * LIMB_BASE + b;
	*low = (digitpair)(hilo & DIGIT_PAIR_MASK);
	return (limb_t)(hilo >> DIGIT_PAIR_SHIFT);
}

/* given a nonnegative integer n representing a number of PyLong digits,
   return an upper bound for the number of limbs required to hold that many
   digits.  Return PY_SSIZE_T_MAX on overflow. */

extern Py_ssize_t
limbsize_from_longsize(Py_ssize_t n) {
	/* compute ceiling(146*n/291); 146/291 is a touch larger than
	   log(2**15)/log(10**9).  we're making n smaller, so there's no
	   danger of overflow. */
	if (n < 10000000)
		/* n not too large (usual case) */
		return (n*146+291 - 1)/291;
	else
		/* n large; rewrite the above expression to avoid overflow */
		return n/291*146 + (n%291*146 + 290)/291;
}

extern Py_ssize_t
longsize_from_limbsize(Py_ssize_t n) {
	/* here we need a rational upper bound for log(10^9)/log(2^30) = 
	   0.9965784284...
	   1 would work;  874/877 =
	   0.9965792474...
	   is slightly better.  Let's go with 1 and waste some space... */
	return n;
}

extern unsigned long
limb_hash(limb_t x) {
	return (unsigned long)x;
}
