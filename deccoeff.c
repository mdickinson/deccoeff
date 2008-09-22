/*
 * Module implementing arbitrary-precision natural number arithmetic
 * in a decimal base.  As the name 'deccoeff' suggests, these numbers
 * are intended to be used as the coefficients for Decimal instances.
 *
 * Author: Mark Dickinson.  Licensed to the PSF under a Contributor
 * Agreement.
 */

/*
 *  To do
 *  -----
 *  write Deccoeff-specific tests
 *  improve and correct documentation
 *  add configuration to detect presence of stdint.h and stdbool.h
 *    and make appropriate substitutions if necessary
 *  fast recursive algorithms for multiplication, division, base conversion
 *  fix usage of stdint and stdbool
 *  implement two and three-argument pow
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Python.h"
#include "longintrepr.h"
#include "string.h"
#include "config.h"

/* substitute for bool if necessary */

#ifdef HAVE_STDBOOL_H
#  include <stdbool.h>
#else
#  ifndef HAVE__BOOL
#    define _Bool int;
#  endif
#  define bool _Bool;
#  define false 0;
#  define true 1;
#endif

/* limb_t is an integer type that can hold any integer in the range [0,
   LIMB_BASE).  double_limb_t should be able to hold any integer in the range
   [0, LIMB_BASE*LIMB_BASE), while digit_limb_t should be able to hold any
   integer in the range [0, LIMB_BASE*PyLong_BASE).

   BASEC_P/BASEC_Q is an upper bound for log(PyLong_BASE)/log(LIMB_BASE);
   BASECI_P/BASECI_Q is an upper bound for log(LIMB_BASE)/log(PyLong_BASE). */

/* Decide which types to use for limb_t and double_limb_t */

#if defined(HAVE_STDINT_H)
#include <stdint.h>
#endif

#if defined(HAVE_INTTYPES_H)
#include <inttypes.h>
#endif

#if (defined(INT32_MAX) || defined(int32_t)) && \
	(defined(INT64_MAX) || defined(int64_t))

/* use int32_t and int64_t if available, with 9 digits to a limb... */

typedef int32_t limb_t;
typedef int64_t double_limb_t;
typedef int64_t digit_limb_t;
#define LIMB_DIGITS 9
#define LIMB_MAX ((limb_t)999999999)
#define BASEC_P 5553
#define BASEC_Q 11068
#define BASECI_P 4369
#define BASECI_Q 2192
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

#else

/* ... else fall back to short and long, and make LIMB_DIGITS 4. */

typedef short limb_t;
typedef long double_limb_t;
typedef long digit_limb_t;
#define LIMB_DIGITS 4
#define LIMB_MAX ((limb_t)9999)
#define BASEC_P 5890
#define BASEC_Q 6649
#define BASECI_P 1717
#define BASECI_Q 1521
static limb_t powers_of_ten[LIMB_DIGITS+1] = {
	1,
	10,
	100,
	1000,
	10000
};

#endif

#define LIMB_ZERO ((limb_t)0)
#define LIMB_ONE ((limb_t)1)
#define LIMB_BASE (LIMB_MAX+1)


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

/* XXX replace limb_error with Py_FataError eventually */

void
limb_error(const char *msg)
{
	fprintf(stderr, "%s\n", msg);
	abort();
}

/* add with carry: compute a + b + c, put result in *r and return
   the new carry. */

static bool
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

static bool
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

static bool
limb_eq(limb_t a, limb_t b)
{
	return a == b;
}

/* multiply and add with two addends; this is useful for long
   multiplication.  Store the low part of the result in *low, and return
   the high part. */

static limb_t
limb_fmaa(limb_t *low, limb_t a, limb_t b, limb_t c, limb_t d) {
	double_limb_t hilo;
	hilo = (double_limb_t)a * b + c + d;
	*low = (limb_t)(hilo%LIMB_BASE);
	return (limb_t)(hilo/LIMB_BASE);
}

/* division: divide high*BASE+low by c, giving a quotient (returned)
   and a remainder (stored in *rem).  Requires that high < c (and
   hence that c is nonzero). */

static limb_t
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

static Py_ssize_t
limb_dsr(limb_t x) {
	Py_ssize_t i;
	if (x == 0)
		limb_error("0 has no most significant digit");
	for (i=0; powers_of_ten[i] <= x; i++);
	return i;
}

/* shift a limb right n places (0 <= n <= LIMB_DIGITS) */

static limb_t
limb_sar(limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_sar");
	return x / powers_of_ten[n];
}

/* rotate right and split; like limb_sar, but also puts piece that was shifted
   out into the top of *res */

static limb_t
limb_split(limb_t *res, limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_split");
	*res = x % powers_of_ten[n] * powers_of_ten[LIMB_DIGITS-n];
	return x / powers_of_ten[n];
}

/* rotate left and split */

static limb_t
limb_splitl(limb_t *res, limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_splitl");
	*res = x / powers_of_ten[LIMB_DIGITS-n];
	return x % powers_of_ten[LIMB_DIGITS-n] * powers_of_ten[n];
}

/* shift a limb left n places (0 <= n <= LIMB_DIGITS) */

static limb_t
limb_sal(limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_sal");
	return x % powers_of_ten[LIMB_DIGITS-n] * powers_of_ten[n];
}

/* select the bottom n digits of a limb (0 <= n <= LIMB_DIGITS) */

static limb_t
limb_low(limb_t x, Py_ssize_t n) {
	if (!(0 <= n && n <= LIMB_DIGITS))
		limb_error("invalid shift count in limb_low");
	return x % powers_of_ten[n];
}

static limb_t
limb_shift_digit_in(limb_t x, int i) {
	if (!(0 <= i && i <= 9))
		limb_error("invalid digit in limb_shift_digit_in");
	if (x >= LIMB_BASE/10)
		x %= LIMB_BASE/10;
	return x * 10 + (limb_t)i;
}

static int
limb_shift_digit_out(limb_t *x, limb_t a) {
	int result;
	result = a % 10;
	*x = a / 10;
	return result;
}

/* type digit_limb_t is an integer type large enough to hold
   any value in the range [0, LIMB_BASE * PyLong_BASE] */

/* given limb a and digit pair b, write a*2**30 + b in the form
   c*LIMB_BASE + d, with c a digit pair and d a limb.  Return c and
   put d in *low. */

static digit
limb_digit_swap(limb_t *low, limb_t a, digit b)
{
	digit_limb_t hilo;
	hilo = ((digit_limb_t)a << PyLong_SHIFT) + b;
	*low = (limb_t)(hilo%LIMB_BASE);
	return (digit)(hilo/LIMB_BASE);
}

/* reverse of the above: given a digit pair a and limb b, write LIMB_BASE*a +
   b in the form c*(2**30) + d, c a limb and d a digit pair.  Return c and put
   d in *low. */

static limb_t
digit_limb_swap(digit *low, digit a, limb_t b)
{
	digit_limb_t hilo;
	hilo = (digit_limb_t)a * LIMB_BASE + b;
	*low = (digit)(hilo & PyLong_MASK);
	return (limb_t)(hilo >> PyLong_SHIFT);
}

static unsigned long
limb_hash(limb_t x) {
	return (unsigned long)x;
}

/* derived operations */

/* a < b iff a - b overflows */

static bool
limb_lt(limb_t a, limb_t b)
{
	limb_t dummy;
	return limb_sbb(&dummy, a, b, false);
}

/* a <= b iff a - b - 1 overflows */

static bool
limb_le(limb_t a, limb_t b)
{
	limb_t dummy;
	return limb_sbb(&dummy, a, b, true);
}

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

/*
   deccoeffs are represented in base BASE, for BASE a suitable power of 10.
   I'll use the word 'limb' to refer to a base BASE digit, and 'digit' to
   refer to a usual decimal digit in the range 0 through 9.  Thus a single
   limb holds some constant number of digits; that constant is called
   LIMB_DIGITS below.

   The C type 'limb' should be large enough to hold a single digit in
   this base: i.e. all numbers in the range [0, BASE).  Type
   'double_limb' should be large enough to hold all numbers in the
   range [0, BASE*BASE).
*/

#define MODULE_NAME "_decimal"
#define CLASS_NAME "Deccoeff"

/*********************************
 * Arithmetic on arrays of limbs *
 *********************************/

/* Low-level functions for operating on arrays of limbs.  These functions
   don't take care of memory allocation; they assume that sufficient space is
   provided for their results. */

/* increment n-limb number a if carry is true, else just copy it; gives n-limb
   result and returns a carry */

static bool
limbs_incc(limb_t *res, const limb_t *a, Py_ssize_t a_size, bool carry)
{
	Py_ssize_t i;
	for (i=0; i < a_size; i++)
		carry = limb_adc(res+i, a[i], LIMB_ZERO, carry);
	return carry;
}

/* decrement n-limb number a if carry is true, else just copy it; gives n-limb
   result and returns a carry */

static bool
limbs_decc(limb_t *res, const limb_t *a, Py_ssize_t a_size, bool carry)
{
	Py_ssize_t i;
	for (i=0; i < a_size; i++)
		carry = limb_sbb(res+i, a[i], LIMB_ZERO, carry);
	return carry;
}

/* add n-limb numbers a and b, producing an n-limb result res and a carry */

static bool
limbs_add(limb_t *res, const limb_t *a, const limb_t *b, Py_ssize_t n)
{
	Py_ssize_t i;
	bool carry;
	carry = false;
	for (i=0; i < n; i++)
		carry = limb_adc(res+i, a[i], b[i], carry);
	return carry;
}

/* subtract n-limb numbers a and b, giving n-limb difference res and a
   carry */

static bool
limbs_sub(limb_t *res, const limb_t *a, const limb_t *b, Py_ssize_t n)
{
	Py_ssize_t i;
	bool carry;
	carry = false;
	for (i=0; i < n; i++)
		carry = limb_sbb(res+i, a[i], b[i], carry);
	return carry;
}

/* multiply m-limb number a by single limb x, getting m-limb result res and
   an extra high limb (returned separately). */

static limb_t
limbs_mul1(limb_t *res, const limb_t *a, Py_ssize_t a_size, limb_t x)
{
	limb_t high;
	Py_ssize_t i;
	high = LIMB_ZERO;
	for (i=0; i < a_size; i++)
		high = limb_fmaa(res+i, a[i], x, high, LIMB_ZERO);
	return high;
}

/* multiply m-limb a by n-limb b, getting m+n-limb result res */

static void
limbs_mul(limb_t *res, const limb_t *a, Py_ssize_t a_size,
	  const limb_t *b, Py_ssize_t b_size)
{
	Py_ssize_t i, j;
	limb_t hiword;
	for (j=0; j < b_size; j++)
		res[j] = LIMB_ZERO;
	for (i=0; i < a_size; i++) {
		hiword = LIMB_ZERO;
		for (j=0; j < b_size; j++)
			hiword = limb_fmaa(res+i+j, a[i], b[j],
					   res[i+j], hiword);
		res[i+j] = hiword;
	}
}

/* divide m-limb number a by single limb x, giving m-limb quotient res
   and returning the (single limb) remainder */

static limb_t
limbs_div1(limb_t *res, const limb_t *a, Py_ssize_t m, limb_t x)
{
	limb_t high;
	Py_ssize_t i;
	high = LIMB_ZERO;
	for (i = m-1; i >= 0; i--)
		res[i] = limb_div(&high, high, a[i], x);
	return high;
}

/* divide m-limb a by n-limb b, giving an (m-n+1)-limb quotient and
   n-limb remainder.  Assumes that the top limb of b is nonzero and
   that m >= n.  w provides m+n+1 limbs of workspace. */

static void
limbs_div(limb_t *quot, limb_t *rem, const limb_t *a, Py_ssize_t m, const limb_t *b,
	  Py_ssize_t n, limb_t *w)
{
	limb_t scale, top, a_top, b_top, q, dummy;
	limb_t *aa, *bb;
	bool carry;
	Py_ssize_t j;

	/* top limb of b should be nonzero; a should contain at least as many
	   limbs as b */
	assert(m >= n && !limb_eq(b[n-1], LIMB_ZERO));

	/* compute scale factor for normalization: floor(LIMB_BASE /
	   (b_top+1)) */
	carry = limb_adc(&scale, b[n-1], LIMB_ONE, false);
	if (carry)
		scale = LIMB_ONE;
	else
		scale = limb_div(&scale, LIMB_ONE, LIMB_ZERO, scale);

	/* scale a and b */
	top = limbs_mul1(w, b, n, scale);
	bb = w;
	assert(limb_eq(top, LIMB_ZERO));

	top = limbs_mul1(w+n, a, m, scale);
	aa = w+n;

	/* catch most cases where quotient only needs m-n limbs */
	if (limb_eq(top, LIMB_ZERO) && limb_lt(aa[m-1], bb[n-1]))
		quot[m-n] = LIMB_ZERO;
	else {
		aa[m] = top;
		m++;
	}

	b_top = bb[n-1];
	aa += m-n;
	for (j = m-n-1; j >= 0; j--) {
		aa--;
		a_top = aa[n];
		assert(limb_le(a_top, b_top));
		/* quotient q = aa / bb; may be overestimate */
		if (limb_eq(a_top, b_top))
			q = LIMB_MAX;
		else
			q = limb_div(&dummy, a_top, aa[n-1], b_top);
		/* compute bottom n limbs of aa[j:] - q*bb */
		top = limbs_mul1(rem, bb, n, q);
		carry = limbs_sub(aa, aa, rem, n);
		carry = limb_adc(&top, top, LIMB_ZERO, carry);
		assert(!carry);
		assert(limb_le(a_top, top));
		/* correct if necessary */
		while (limb_lt(a_top, top)) {
			carry = limbs_add(aa, aa, bb, n);
			carry = limb_adc(&a_top, a_top, LIMB_ZERO, carry);
			assert(!carry);
			carry = limb_sbb(&q, q, LIMB_ONE, false);
			assert(!carry);
		}
		quot[j] = q;
	}
	/* OPT: optimization opportunity---this final step
	   can be skipped when all we want is the quotient. */
	top = limbs_div1(rem, aa, n, scale);
	assert(limb_eq(top, LIMB_ZERO));
}

/* shift m-limb number right n digits (not limbs!).  i.e., divide by
   10**n.  n should satisfy 0 <= n <= LIMB_DIGITS*m, and the result
   res has m - n/LIMB_DIGITS limbs. */

static void
limbs_rshift(limb_t *res, const limb_t *a, Py_ssize_t m, Py_ssize_t n)
{
	Py_ssize_t n_limbs, n_digits, i;
	limb_t limb_top, limb_bot, rem;
	bool carry;

	assert(n >= 0);
	assert(n <= m*LIMB_DIGITS);

	n_limbs = n / LIMB_DIGITS;
	n_digits = n % LIMB_DIGITS;
	if (n_digits == 0) {
		for (i = 0; i < m - n_limbs; i++)
			res[i] = a[i+n_limbs];
		return;
	}
	rem = limb_sar(a[n_limbs], n_digits);
	for (i = 0; i < m-n_limbs-1; i++) {
		limb_top = limb_split(&limb_bot, a[n_limbs+i+1], n_digits);
		carry = limb_adc(res+i, rem, limb_bot, false);
		assert(!carry);
		rem = limb_top;
	}
	res[i] = rem;
}

/* get slice a[m:n] of (the decimal digits of) an integer a, from
   digit m up to (but not including) digit n, giving a result res with
   1+(n-m-1)/LIMB_DIGITS limbs.  Assumes that 0 <= m < n and that a
   has at least 1+(n-1)/LIMB_DIGITS limbs. */

static void
limbs_slice(limb_t *res, const limb_t *a, Py_ssize_t m, Py_ssize_t n)
{
	Py_ssize_t mlimbs, mdigits, res_limbs, res_digits, i;
	limb_t out, limb_bot, limb_top;
	bool carry;
	mlimbs = m / LIMB_DIGITS;
	mdigits = m % LIMB_DIGITS;
	res_limbs = (n-1-m) / LIMB_DIGITS;
	res_digits = (n-1-m) % LIMB_DIGITS + 1;
	out = limb_sar(a[mlimbs++], mdigits);
	for (i = 0; i < res_limbs; i++) {
		limb_top = limb_split(&limb_bot, a[mlimbs++], mdigits);
		carry = limb_adc(res+i, out, limb_bot, false);
		assert(!carry);
		out = limb_top;
	}
	if (res_digits > LIMB_DIGITS - mdigits) {
		limb_bot = limb_sal(a[mlimbs++], LIMB_DIGITS - mdigits);
		carry = limb_adc(&out, out, limb_bot, false);
		assert(!carry);
	}
	res[i] = limb_low(out, res_digits);
}

/* shift m-limb number left n digits (shifting zeros in); i.e.,
   multiply by 10**n.  Assumes res has size m +
   1+(n-1)/LIMB_DIGITS (or just size m if n == 0). */

static void
limbs_lshift(limb_t *res, const limb_t *a, Py_ssize_t m, Py_ssize_t n)
{
	Py_ssize_t n_limbs, n_digits, i;
	limb_t limb_top, tmp;
	bool carry;

	assert(n >= 0);
	n_limbs = n / LIMB_DIGITS;
	n_digits = n % LIMB_DIGITS;
	for (i = 0; i < n_limbs; i++)
		res[i] = LIMB_ZERO;
	if (n_digits == 0) {
		for (; i < n_limbs+m; i++)
			res[i] = a[i-n_limbs];
		return;
	}
	res[i++] = limb_splitl(&limb_top, a[0], n_digits);
	for (; i < n_limbs+m; i++) {
		carry = limb_adc(res+i, limb_top,
				 limb_splitl(&tmp, a[i-n_limbs], n_digits),
				 false);
		assert(!carry);
		limb_top = tmp;
	}
	res[i] = limb_top;
}

/* Conversion to and from strings */

/* Convert a character array to an array of limbs.  Returns false on success,
   true if any of the characters was not a valid digit; in the latter case,
   the contents of a are undefined. a should provide ceiling(slen /
   LIMB_DIGITS) limbs. */

static bool
limbs_from_string(limb_t *a, const char *s, Py_ssize_t s_len)
{
	Py_ssize_t i, k, digits_in_limb;
	limb_t acc;
	char c;

	k = (s_len+LIMB_DIGITS-1) / LIMB_DIGITS;
	digits_in_limb = (s_len+LIMB_DIGITS-1) % LIMB_DIGITS + 1;
	acc = LIMB_ZERO;
	for (i = 0; i < s_len; i++) {
		c = s[i];
		if (c < '0' || c > '9')
			return true;
		acc = limb_shift_digit_in(acc, (int)(c - '0'));
		digits_in_limb--;
		if (digits_in_limb == 0) {
			digits_in_limb = LIMB_DIGITS;
			a[--k] = acc;
			acc = LIMB_ZERO;
		}
	}
	assert(digits_in_limb = LIMB_DIGITS);
	return false;
}

/* Base conversion, from base 2**15 to base LIMB_BASE.

   Convert an array of (base 2**15) digits for a Python long to an
   array of limbs representing the same number.  Returns the number of
   limbs of a filled.  The result is normalized, in the sense that if
   the returned size a_size is nonzero then a[a_size-1] is nonzero.

   a should have at least:

     ceiling(b_size * log(2**15)/log(LIMB_BASE))

   limbs available.
 */

static Py_ssize_t
limbs_from_longdigits(limb_t *a, const digit *b, Py_ssize_t b_size)
{
	Py_ssize_t i, j, a_size;
	digit high;
	a_size = 0;
	for (j = b_size-1; j >= 0; j--) {
		high = b[j];
		for (i=0; i < a_size; i++)
			high = limb_digit_swap(a+i, a[i], high);
		while (high != 0) {
			high = limb_digit_swap(a+a_size, LIMB_ZERO, high);
			a_size++;
		}
	}
	return a_size;
}

/* base conversion, from base LIMB_BASE to base 2**30. */

static Py_ssize_t
limbs_to_longdigits(digit *b, const limb_t *a, Py_ssize_t a_size)
{
	Py_ssize_t i, j, b_size;
	limb_t high;
	b_size = 0;
	for (i = a_size-1; i >= 0; i--) {
		high = a[i];
		for (j = 0; j < b_size; j++)
			high = digit_limb_swap(b+j, b[j], high);
		while (!limb_eq(high, LIMB_ZERO)) {
			high = digit_limb_swap(b+b_size, 0, high);
			b_size++;
		}
	}
	return b_size;
}

/**************************
 * deccoeff : definitions *
 **************************/

/* We place an upper bound MAX_DIGITS on the number of decimal digits (*not*
   the number of limbs) in a deccoeff.  MAX_DIGITS should fit into a
   Py_ssize_t, so that length and indexing always make sense.  Here we take
   MAX_DIGITS = 10**9; then assuming that a Py_ssize_t is at least 32 bits
   we're safe from overflow even when adding two digit counts. */

#define MAX_DIGITS 1000000000

/* N.B. We really want ob_limbs[0] in this definition, but that isn't
   standards-compliant C (though gcc permits it).  This makes
   computing tp_basicsize a little tricky. */

typedef struct {
  PyObject_VAR_HEAD
  limb_t ob_limbs[1];
} deccoeff;

/* XXX this results in wasted memory (DECCOEFF_BASICSIZE is too
   large); fix me! */

#define DECCOEFF_ITEMSIZE sizeof(limb_t)
#define DECCOEFF_BASICSIZE sizeof(deccoeff)

static PyTypeObject deccoeff_DeccoeffType;

/* allocate a new decimal integer with 'size' limbs */

static deccoeff *
_deccoeff_new(Py_ssize_t size)
{
	return PyObject_NEW_VAR(deccoeff, &deccoeff_DeccoeffType, size);
}

/* normalize a deccoeff, by setting the size correctly. */

static deccoeff *
deccoeff_normalize(deccoeff *v)
{
	Py_ssize_t v_size;
	v_size = Py_SIZE(v);
	while (v_size > 0 && limb_eq(v->ob_limbs[v_size-1], LIMB_ZERO))
		--v_size;
	Py_SIZE(v) = v_size;
	return v;
}

/* returns its argument (with reference count unaffected) if the argument has
   MAX_DIGITS or fewer digits.  Otherwise it decrements the reference count
   for its argument, sets OverflowError, and returns NULL. */

static deccoeff *
deccoeff_checksize(deccoeff *v)
{
	Py_ssize_t v_size, topdigits;
	bool small;
	v_size = Py_SIZE(v);

	/* something with MAX_DIGITS digits has exactly
	   (MAX_DIGITS-1)/LIMB_DIGITS+1 limbs;  its top limb has
	   (MAX_DIGITS-1)%LIMB_DIGITS+1 digits */
	if (v_size < (MAX_DIGITS-1)/LIMB_DIGITS+1)
		small = true;
	else if (v_size == (MAX_DIGITS-1)/LIMB_DIGITS+1) {
		topdigits = limb_dsr(v->ob_limbs[v_size-1]);
		small = topdigits <= (MAX_DIGITS-1)%LIMB_DIGITS+1;
	}
	else
		small = false;

	if (small)
		return v;
	Py_DECREF(v);
	PyErr_SetString(PyExc_OverflowError,
			"Deccoeff instance has too many digits");
	return NULL;
}

/* return a new reference to the zero deccoeff */

static deccoeff *
deccoeff_zero(void)
{
	return _deccoeff_new(0);
}

static deccoeff *
deccoeff_from_string_and_size(const char *s, Py_ssize_t s_len) {
	Py_ssize_t z_size;
	deccoeff *z;
	bool invalid;

	if (s_len > MAX_DIGITS) {
		PyErr_SetString(PyExc_OverflowError,
				"too many digits");
		return NULL;
	}

	z_size = (s_len + LIMB_DIGITS - 1) / LIMB_DIGITS;
	z = _deccoeff_new(z_size);
	if (z == NULL)
		return NULL;

	invalid = limbs_from_string(z->ob_limbs, s, s_len);
	if (invalid) {
		Py_DECREF(z);
		PyErr_SetString(PyExc_ValueError,
				"nondigit character in input");
		return NULL;
	}
	return deccoeff_normalize(z);
}




/***************************
 * Arithmetic on deccoeffs *
 ***************************/

/* General rules: if the result of any arithmetic operation falls
   outside the range [0, 10**MAX_DIGITS) then OverflowError is raised.
   Results are always normalized. */


/* convert an arbitrary PyObject to a deccoeff.  If conversion is implemented
   for this type, return false and put the result of the attempted conversion
   (which may be NULL) in *z.   Otherwise, return true. */

static deccoeff *deccoeff_from_PyLong(PyLongObject *a);

static bool
convert_to_deccoeff(deccoeff **z, PyObject *v)
{
	PyLongObject *w;
	if (v->ob_type == &deccoeff_DeccoeffType) {
		Py_INCREF(v);
		*z = (deccoeff *)v;
		return false;
	}
	else if (PyIndex_Check(v)) {
		w = (PyLongObject *)PyNumber_Index(v);
		if (w == NULL) {
			*z = NULL;
			return false;
		}
		*z = deccoeff_from_PyLong(w);
		Py_DECREF(w);
		return false;
	}
	else
		return true;
}


#define DECCOEFF_WRAP_BINOP(PO_func, DC_func)   \
static PyObject *                               \
PO_func(PyObject *v, PyObject *w)               \
{                                               \
	deccoeff *a, *b;                        \
	PyObject *z;                            \
	if (convert_to_deccoeff(&a, v)) {       \
		Py_INCREF(Py_NotImplemented);   \
		return Py_NotImplemented;       \
	}                                       \
	if (convert_to_deccoeff(&b, w)) {       \
		Py_DECREF(a);                   \
		Py_INCREF(Py_NotImplemented);   \
		return Py_NotImplemented;       \
	}                                       \
	z = (PyObject *)(DC_func(a, b));        \
	Py_DECREF(a);                           \
	Py_DECREF(b);                           \
	return z;                               \
}

/* addition */

static deccoeff *
_deccoeff_add(deccoeff *a, deccoeff *b)
{
	Py_ssize_t a_size, b_size;
	deccoeff *z;
	bool carry;
	a_size = Py_SIZE(a);
	b_size = Py_SIZE(b);
	if (a_size < b_size) {
		deccoeff *temp;
		Py_ssize_t temp_size;
		temp = a; a = b; b = temp;
		temp_size = a_size; a_size = b_size; b_size = temp_size;
	}
	z = _deccoeff_new(a_size+1);
	if (z == NULL)
		return z;
	carry = limbs_add(z->ob_limbs, a->ob_limbs, b->ob_limbs, b_size);
	carry = limbs_incc(z->ob_limbs+b_size, a->ob_limbs+b_size,
		   a_size-b_size, carry);
	z->ob_limbs[a_size] = carry ? LIMB_ONE : LIMB_ZERO;
	return deccoeff_checksize(deccoeff_normalize(z));
}

/* subtraction */

static deccoeff *
_deccoeff_subtract(deccoeff *a, deccoeff *b)
{
	Py_ssize_t a_size, b_size;
	deccoeff *z;
	bool carry;
	a_size = Py_SIZE(a);
	b_size = Py_SIZE(b);
	if (a_size < b_size)
		goto Overflow;
	z = _deccoeff_new(a_size);
	if (z==NULL)
		return NULL;
	carry = limbs_sub(z->ob_limbs, a->ob_limbs, b->ob_limbs, b_size);
	carry = limbs_decc(z->ob_limbs+b_size, a->ob_limbs+b_size,
			   a_size-b_size, carry);
	if (!carry)
		return deccoeff_normalize(z);
	Py_DECREF(z);
  Overflow:
	PyErr_SetString(PyExc_OverflowError, "difference is negative");
	return NULL;
}

/* multiplication */

static deccoeff *
_deccoeff_multiply(deccoeff *a, deccoeff *b)
{
	Py_ssize_t a_size, b_size;
	deccoeff *z;
	a_size = Py_SIZE(a);
	b_size = Py_SIZE(b);
	z = _deccoeff_new(a_size + b_size);
	if (z == NULL)
		return NULL;
	limbs_mul(z->ob_limbs, a->ob_limbs, a_size, b->ob_limbs, b_size);
	return deccoeff_checksize(deccoeff_normalize(z));
}

/* division of a by b: returns quotient and puts remainder in *r.  On
   failure, both the quotient and the remainder are NULL.  Raises
   ZeroDivisionError on division by zero. */

static deccoeff *
_deccoeff_division(deccoeff **r, deccoeff *a, deccoeff *b) {
	deccoeff *w, *rem, *quot;
	Py_ssize_t a_size, b_size;
	a_size = Py_SIZE(a);
	b_size = Py_SIZE(b);
	/* special cases b == 0, a_size < b_size */
	if (b_size == 0) {
		PyErr_SetString(PyExc_ZeroDivisionError,
				"division by zero");
		*r = NULL;
		return NULL;
	}
	if (a_size < b_size) {
		Py_INCREF(a);
		*r = a;
		return deccoeff_zero();
	}
	quot = _deccoeff_new(a_size + 1 - b_size);
	if (quot == NULL) {
		*r = NULL;
		return NULL;
	}
	rem = _deccoeff_new(b_size);
	if (rem == NULL) {
		Py_DECREF(quot);
		*r = NULL;
		return NULL;
	}
	if (b_size == 1)
		/* fast path for division by a single limb */
		rem->ob_limbs[0] = limbs_div1(quot->ob_limbs, a->ob_limbs,
					      a_size, b->ob_limbs[0]);
	else {
		/* long division */
		w = _deccoeff_new(a_size + 1 + b_size);
		if (w == NULL) {
			Py_DECREF(rem);
			Py_DECREF(quot);
			*r = NULL;
			return NULL;
		}
		limbs_div(quot->ob_limbs, rem->ob_limbs, a->ob_limbs, a_size,
			  b->ob_limbs, b_size, w->ob_limbs);
		Py_DECREF(w);
	}
	*r = deccoeff_normalize(rem);
	return deccoeff_normalize(quot);
}

/* remainder: raises ZeroDivisionError if b is zero */

static deccoeff *
_deccoeff_remainder(deccoeff *a, deccoeff *b) {
	deccoeff *quot, *rem;
	quot = _deccoeff_division(&rem, a, b);
	if (rem == NULL)
		return NULL;
	Py_DECREF(quot);
	return rem;
}

/* divmod: raises ZeroDivisionError if b is zero */

static PyObject *
_deccoeff_divmod(deccoeff *a, deccoeff *b) {
	deccoeff *quot, *rem;

	quot = _deccoeff_division(&rem, a, b);
	if (quot == NULL)
		return NULL;
	return Py_BuildValue("OO", (PyObject *)quot, (PyObject *)rem);
}

/* negation: succeeds only for a == 0; for anything else, OverflowError is
   raised */

static deccoeff *
deccoeff_negative(deccoeff *a)
{
	Py_ssize_t a_size;

	a_size = Py_SIZE(a);
	if (a_size == 0)
		return deccoeff_zero();
	PyErr_SetString(PyExc_OverflowError, "negation is negative");
	return NULL;
}

/* unary plus: does nothing */

static deccoeff *
deccoeff_positive(deccoeff *a)
{
	Py_INCREF(a);
	return a;
}

/* bool: return True if a is nonzero, else False */

static int
deccoeff_bool(deccoeff *a)
{
	return Py_SIZE(a) != 0;
}

/* left shift; second operand is a Python integer, not a
   deccoeff. raises ValueError if second operand is negative */

/* how should this behave with respect to other types?
   i.e., should 2 << DI('3') be an error?
   For now, yes:  can't see much value in 2 << DI('3').

   So: first argument should be a Deccoeff.  Second argument
   may be a decoeff, may be something else that can be interpreted as
   an index.
 */

static PyObject *
deccoeff_lshift(PyObject *v, PyObject *b) {
	Py_ssize_t n, a_size;
	deccoeff *z, *a;

	/* shifting another type by a deccoeff is not supported */
	if (v->ob_type != &deccoeff_DeccoeffType || !PyIndex_Check(b)) {
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	}
	a = (deccoeff *)v;

	/* attempt to interpret b as a nonnegative index */
	n = PyNumber_AsSsize_t(b, NULL);
	if (n == -1 && PyErr_Occurred())
		return NULL;
	else if (n < 0) {
		PyErr_SetString(PyExc_ValueError, "negative shift count");
		return NULL;
	}

	a_size = Py_SIZE(a);
	if (a_size == 0)
		return (PyObject *)deccoeff_zero();
	if (n >= MAX_DIGITS) {
		PyErr_SetString(PyExc_OverflowError,
				"Deccoeff instance has too many digits");
		return NULL;
	}
	z = _deccoeff_new(a_size + (n+LIMB_DIGITS-1) / LIMB_DIGITS);
	if (z==NULL)
		return NULL;
	limbs_lshift(z->ob_limbs, a->ob_limbs, a_size, n);
	return (PyObject *)deccoeff_checksize(deccoeff_normalize(z));
}

/* right shift; second operand is a Python integer, not a deccoeff.
   Raises ValueError if second operand is negative. */

static PyObject *
deccoeff_rshift(PyObject *v, PyObject *b) {
	Py_ssize_t n, a_size, shift;
	deccoeff *z, *a;

	/* shifting another type by a deccoeff is not supported */
	if (v->ob_type != &deccoeff_DeccoeffType || !PyIndex_Check(b)) {
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	}
	a = (deccoeff *)v;

	/* attempt to interpret b as a nonnegative index */
	n = PyNumber_AsSsize_t(b, NULL);
	if (n == -1 && PyErr_Occurred())
		return NULL;
	else if (n < 0) {
		PyErr_SetString(PyExc_ValueError, "negative shift count");
		return NULL;
	}

	a_size = Py_SIZE(a);
	shift = n / LIMB_DIGITS;
	if (shift >= a_size)
		return (PyObject *)deccoeff_zero();
	z = _deccoeff_new(a_size - shift);
	if (z==NULL)
		return NULL;
	limbs_rshift(z->ob_limbs, a->ob_limbs, a_size, n);
	return (PyObject *)deccoeff_normalize(z);
}

/* slice: the slice indices should be nonnegative integers;  the step
   argument is not supported and is expected to be None. */

static deccoeff *
deccoeff_subscript(deccoeff *a, PyObject *b)
{
	PySliceObject *slice;
	Py_ssize_t a_size, start, stop, defstop;
	deccoeff *z;

	a_size = Py_SIZE(a);
	defstop = a_size * LIMB_DIGITS;
	if (PySlice_Check(b)) {
		/* for a slice, extract start, stop and step.  We'd like to
		use PySlice_GetIndicesEx here, but we can't supply a sensible
		length.  So we imitate it instead, which is a little bit
		naughty since it involves using the private function
		_PyEval_SliceIndex. */
		slice = (PySliceObject *)b;
		if (slice->step != Py_None) {
			PyErr_SetString(PyExc_ValueError,
					"step unsupported in deccoeff slice");
			return NULL;
		}
		if (slice->start == Py_None)
			start = 0;
		else {
			if (!_PyEval_SliceIndex(slice->start, &start))
				return NULL;
			if (start < 0) {
				PyErr_SetString(PyExc_ValueError,
						"start should be nonnegative");
				return NULL;
			}
		}
		defstop = a_size * LIMB_DIGITS;
		if (slice->stop == Py_None)
			stop = defstop;
		else {
			if (!_PyEval_SliceIndex(slice->stop, &stop))
				return NULL;
			if (stop < 0) {
				PyErr_SetString(PyExc_ValueError,
						"stop should be nonnegative");
				return NULL;
			}
		}
	}
	else {
		start = PyNumber_AsSsize_t(b, NULL);
		if (start == -1 && PyErr_Occurred())
			return NULL;
		else if (start < 0) {
			PyErr_SetString(PyExc_ValueError,
					"index should be nonnegative");
			return NULL;
		}
		/* clip large results to MAX_DIGITS */
		if (start > MAX_DIGITS)
			start = MAX_DIGITS;
		stop = start + 1;
	}
	/* reduce to case 0 <= start < stop <= a_size * LIMB_DIGITS */
	if (stop > defstop)
		stop = defstop;
	if (stop <= start)
		return deccoeff_zero();

	/* allocate space for result */
	z = _deccoeff_new((stop-start-1)/LIMB_DIGITS + 1);
	if (z==NULL)
		return NULL;
	limbs_slice(z->ob_limbs, a->ob_limbs, start, stop);
	return deccoeff_normalize(z);
}

/* floor division:  raise ZeroDivisionError if b is 0 */

static deccoeff *
_deccoeff_floor_divide(deccoeff *a, deccoeff *b) {
	deccoeff *quot, *rem;
	quot = _deccoeff_division(&rem, a, b);
	if (quot == NULL)
		return NULL;
	Py_DECREF(rem);
	return quot;
}

/* compare: return -1 if a < b, 0 if a == b, 1 if a > b.  Always
   succeeds. */

static int
_deccoeff_compare(deccoeff *a, deccoeff *b)
{
	Py_ssize_t a_size, b_size, i;
	a_size = Py_SIZE(a);
	b_size = Py_SIZE(b);
	if (a_size != b_size)
		return a_size < b_size ? -1 : 1;
	for (i = a_size-1; i >= 0; i--)
		if (!limb_eq(a->ob_limbs[i], b->ob_limbs[i]))
			return limb_lt(a->ob_limbs[i], b->ob_limbs[i]) ? -1 : 1;
	return 0;
}

/* Python 3.x only understands rich comparisons; deccoeff_richcompare wraps
   _deccoeff_compare. */

static PyObject *
deccoeff_richcompare(PyObject *self, PyObject *other, int op)
{
	PyObject *result;
	result = Py_CmpToRich(op, _deccoeff_compare((deccoeff*)self,
						   (deccoeff*)other));
	return result;
}

/* Compute ceiling(n*p/q) without intermediate overflow.  If the result
   would be larger than PY_SSIZE_T_MAX, return -1.  Assumes that
   n is nonnegative, and that (q-1)*(p+1) <= PY_SSIZE_T_MAX. */

static Py_ssize_t
scale_Py_ssize_t(Py_ssize_t n, int p, int q) {
	Py_ssize_t hi, low;
	assert (n >= 0);
	hi = n/q;
	if (hi > PY_SSIZE_T_MAX/p)
		return -1;
	hi *= p;
	low = (n%q*p+q-1)/q;
	if (hi > PY_SSIZE_T_MAX-low)
		return -1;
	return hi+low;
}

/* Create a deccoeff from a Python integer. */

static deccoeff *
deccoeff_from_PyLong(PyLongObject *a)
{
	Py_ssize_t a_size, z_size;
	deccoeff *z;

	a_size = Py_SIZE(a);
	if (a_size < 0) {
		PyErr_SetString(PyExc_OverflowError,
				"Can't convert negative integer to deccoeff");
		return NULL;
	}

	z_size = scale_Py_ssize_t(a_size, BASEC_P, BASEC_Q);
	if (z_size == -1)
		PyErr_SetString(PyExc_OverflowError,
				"Overflow in int to deccoeff conversion\n");
	z = _deccoeff_new(z_size);
	if (z==NULL)
		return NULL;
	Py_SIZE(z) = limbs_from_longdigits(z->ob_limbs, a->ob_digit, a_size);
	return deccoeff_checksize(z);
}

static PyLongObject *
deccoeff_long(deccoeff *a)
{
	Py_ssize_t a_size, z_size;
	PyLongObject *z;

	a_size = Py_SIZE(a);
	z_size = scale_Py_ssize_t(a_size, BASECI_P, BASECI_Q);
	if (z_size == -1)
		PyErr_SetString(PyExc_OverflowError,
				"Overflow in deccoeff to int conversion\n");
	z = _PyLong_New(z_size);
	if (z == NULL)
		return NULL;
	Py_SIZE(z) = limbs_to_longdigits(z->ob_digit, a->ob_limbs, a_size);
	return z;
}

static PyObject *
deccoeff_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	PyObject *x, *n, *result;
	static char *kwlist[] = {"x", 0};
	char *s;
	Py_ssize_t s_len;

	/* not allowing subtypes */
	assert(type == &deccoeff_DeccoeffType);
	x = NULL;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O:" CLASS_NAME,
					 kwlist, &x))
		return NULL;

	if (x == NULL)
		return (PyObject *)deccoeff_zero();
	else if (PyUnicode_Check(x)) {
		s = _PyUnicode_AsStringAndSize(x, &s_len);
		if (s == NULL)
			return NULL;
		return (PyObject *)deccoeff_from_string_and_size(s, s_len);
	}
	else if (PyIndex_Check(x)) {
		n = PyNumber_Index(x);
		if (n == NULL)
			return NULL;
		result = (PyObject *)deccoeff_from_PyLong((PyLongObject *)n);
		Py_DECREF(n);
		return result;
	}
	else if (PyObject_TypeCheck(x, &deccoeff_DeccoeffType)) {
		/* we're not allowing subtypes, so it should
		   be safe just to return x itself here. */
		Py_INCREF(x);
		return x;
	}
	else {
		PyErr_SetString(PyExc_TypeError,
				"unable to create a deccoeff from this type");
		return NULL;
	}
}


/* number of digits of a deccoeff, or 0 if that deccoeff is zero. */

static Py_ssize_t
deccoeff_length(deccoeff *v)
{
	Py_ssize_t v_size;
	v_size = Py_SIZE(v);
	if (v_size == 0)
		return 0;
	return limb_dsr(v->ob_limbs[v_size-1]) + (v_size-1) * LIMB_DIGITS;
}

static PyObject *
deccoeff_str(deccoeff *v)
{
	Py_ssize_t sz, nlimbs;
	limb_t *limb_pointer, *last_limb, limb_value;
	PyObject *str;
	int i;
	Py_UNICODE *p;

	nlimbs = Py_SIZE(v);
	if (nlimbs == 0) {
		/* return empty string */
		str = PyUnicode_FromUnicode(NULL, 1);
		if (str == NULL)
			return NULL;
		p = PyUnicode_AS_UNICODE(str) + 1;
		*p = '\0';
		*--p = '0';
		return str;
	}

	sz = deccoeff_length(v);

	str = PyUnicode_FromUnicode(NULL, sz);
	if (str == NULL)
		return NULL;
	p = PyUnicode_AS_UNICODE(str) + sz;
	*p = '\0';

	/* fill in digits from right to left;  start with the least
	   significant limb_t */
	limb_pointer = v -> ob_limbs;
	last_limb = limb_pointer + nlimbs - 1;
	while (limb_pointer < last_limb) {
		limb_value = *limb_pointer++;
		for (i=0; i < LIMB_DIGITS; i++)
			*--p = '0' + limb_shift_digit_out(
				&limb_value, limb_value);
	}
	/* most significant limb_t */
	limb_value = *limb_pointer;
	assert(!limb_eq(limb_value, LIMB_ZERO));
	while (!limb_eq(limb_value, LIMB_ZERO)) {
		*--p = '0' + limb_shift_digit_out(
			&limb_value, limb_value);
	}
	return str;
}

static PyObject *
deccoeff_repr(deccoeff *v)
{
	PyObject *strv;
	PyObject *result;
	strv = deccoeff_str(v);
	if (strv == NULL)
		return NULL;
	result = PyUnicode_FromFormat(CLASS_NAME "('%U')", strv);
	Py_DECREF(strv);
	return result;
}

static void
deccoeff_dealloc(PyObject *v)
{
	Py_TYPE(v)->tp_free(v);
}

#define HASH_SHIFT 13
#define HASH_MASK ((1<<HASH_SHIFT) - 1)
#define HASH_START 1887730231
#define HASH_BITS 32

static long
deccoeff_hash(deccoeff *v)
{
	unsigned long x;
	long y;
	limb_t *v_start, *v_top;

	/* We don't bother to make hashes equal for ints and corresponding
	   deccoeffs, since that would slow things down.  In fact, we
	   deliberately make it such that a deccoeff is unlikely to hash to the
	   same value as the corresponding int, to avoid subtle errors.*/

	/* produce a 32-bit hash value */
	v_start = v->ob_limbs;
	v_top = v_start + Py_SIZE(v);
	x = HASH_START;
	while (v_top > v_start) {
		x = ((x & ~HASH_MASK) >> HASH_SHIFT) |
			((x & HASH_MASK) << (HASH_BITS-HASH_SHIFT));
		x ^= limb_hash(*--v_top);
	}

	y = (long)x;
	return (y == -1) ? -2 : y;
}

static PyMethodDef deccoeff_methods[] = {
	{NULL, NULL}
};

static PyMappingMethods deccoeff_as_mapping = {
	(lenfunc)deccoeff_length,               /*mp_length*/
	(binaryfunc)deccoeff_subscript,         /*mp_subscript*/
	0, /*mp_ass_subscript*/
};

DECCOEFF_WRAP_BINOP(deccoeff_add, _deccoeff_add)
DECCOEFF_WRAP_BINOP(deccoeff_subtract, _deccoeff_subtract)
DECCOEFF_WRAP_BINOP(deccoeff_multiply, _deccoeff_multiply)
DECCOEFF_WRAP_BINOP(deccoeff_remainder, _deccoeff_remainder)
DECCOEFF_WRAP_BINOP(deccoeff_divmod, _deccoeff_divmod)
DECCOEFF_WRAP_BINOP(deccoeff_floor_divide, _deccoeff_floor_divide)

static PyNumberMethods deccoeff_as_number = {
	deccoeff_add,                           /*nb_add*/
	deccoeff_subtract,                      /*nb_subtract*/
	deccoeff_multiply,                      /*nb_multiply*/
	deccoeff_remainder,                     /*nb_remainder*/
	deccoeff_divmod,                        /*nb_divmod*/
	0, /*nb_power*/
	(unaryfunc) deccoeff_negative,          /*nb_negative*/
	(unaryfunc) deccoeff_positive,          /*nb_positive*/
	(unaryfunc) deccoeff_positive,          /*nb_absolute*/
	(inquiry) deccoeff_bool,                /*nb_bool*/
	0, /*nb_invert*/
	deccoeff_lshift,                        /*nb_lshift*/
	deccoeff_rshift,                        /*nb_rshift*/
	0, /*nb_and*/
	0, /*nb_xor*/
	0, /*nb_or*/
	(unaryfunc) deccoeff_long,              /*nb_int*/
	(unaryfunc) deccoeff_long,              /*nb_long*/
	0, /*nb_float*/
	0, /*nb_inplace_add*/
	0, /*nb_inplace_subtract*/
	0, /*nb_inplace_multiply*/
	0, /*nb_inplace_remainder*/
	0, /*nb_inplace_power*/
	0, /*nb_inplace_lshift*/
	0, /*nb_inplace_rshift*/
	0, /*nb_inplace_and*/
	0, /*nb_inplace_xor*/
	0, /*nb_inplace_or*/
	deccoeff_floor_divide,                  /*nb_floor_divide*/
	0, /*nb_true_divide*/
	0, /*nb_inplace_floor_divide*/
	0, /*nb_inplace_true_divide*/
	(unaryfunc) deccoeff_long,              /*nb_index*/
};

static PyTypeObject deccoeff_DeccoeffType = {
	PyVarObject_HEAD_INIT(&PyType_Type, 0)
	MODULE_NAME "." CLASS_NAME,             /* tp_name */
	DECCOEFF_BASICSIZE,                     /* tp_basicsize */
	DECCOEFF_ITEMSIZE,                      /* tp_itemsize */
	deccoeff_dealloc,                       /* tp_dealloc */
	0, /* tp_print */
	0, /* tp_getattr */
	0, /* tp_setattr */
	0, /* tp_compare */
	(reprfunc)deccoeff_repr,                /* tp_repr */
	&deccoeff_as_number,                    /* tp_as_number */
	0, /* tp_as_sequence */
	&deccoeff_as_mapping,                   /* tp_as_mapping */
	(hashfunc)deccoeff_hash,                /* tp_hash */
	0, /* tp_call */
	(reprfunc)deccoeff_str,                 /* tp_str */
	0, /* tp_getattro */
	0, /* tp_setattro */
	0, /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT,                     /* tp_flags */
	"Decimal integers",                     /* tp_doc */
	0, /* tp_traverse */
	0, /* tp_clear */
	(richcmpfunc)deccoeff_richcompare,      /* tp_richcompare */
	0, /* tp_weaklistoffset */
	0, /* tp_iter */
	0, /* tp_iternext */
	deccoeff_methods,                       /* tp_methods */
	0, /* tp_members */
	0, /* tp_getset */
	0, /* tp_base */
	0, /* tp_dict */
	0, /* tp_descr_get */
	0, /* tp_descr_set */
	0, /* tp_dictoffset */
	0, /* tp_init */
	0, /* tp_alloc */
	deccoeff_new,                           /* tp_new */
	PyObject_Del,                           /* tp_free */
};

static PyMethodDef deccoeff_module_methods[] = {
	{NULL, NULL}
};

static struct PyModuleDef _decimalmodule = {
	PyModuleDef_HEAD_INIT,
	"_decimal",
	"class for decimal integer arithmetic; support for decimal module",
	-1,
	deccoeff_module_methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC
PyInit__decimal(void)
{
	PyObject *m;
	int check;

	if (PyType_Ready(&deccoeff_DeccoeffType) < 0)
		return NULL;

	m = PyModule_Create(&_decimalmodule);
	if (m == NULL)
		return NULL;

	Py_INCREF(&deccoeff_DeccoeffType);
	check = PyModule_AddObject(m, CLASS_NAME,
				   (PyObject *) &deccoeff_DeccoeffType);
	if (check == -1)
		return NULL;
	check = PyModule_AddIntMacro(m, LIMB_DIGITS);
	if (check == -1)
		return NULL;
	check = PyModule_AddIntMacro(m, MAX_DIGITS);
	if (check == -1)
		return NULL;
	return m;
}