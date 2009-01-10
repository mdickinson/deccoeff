/*

   Add check that MAX_DIGITS/LIMB_DIGITS is small enough:  twice it
   should fit into a Py_ssize_t.

*/

/*
 * This file defines a module deccoeff containing two classes:
 * _Decimal and Deccoeff.
 *
 * deccoeff._Decimal is a skeletal base class for the decimal.Decimal class.
 * As time goes on, the aim is to move more and more code from the Python
 * decimal.py file to the _Decimal class, and eventually rename _Decimal
 * to Decimal.
 *
 * deccoeff.Deccoeff is a class implementing arbitrary-precision unsigned
 * integer arithmetic in a decimal base.  In addition to the usual arithmetic
 * operations, Deccoeff instances support slicing and element access for
 * retrieving individual digits or sequences of digits.
 *
 * As the name suggests, Deccoeff instances are intended to be used as the
 * coefficients for _Decimal instances.
 *
 * Author: Mark Dickinson (dickinsm@gmail.com).
 * Licensed to the PSF under a Contributor Agreement.
 */

#include "Python.h"
#include "longintrepr.h"
#include "deccoeff_config.h"

#include <stddef.h>  /* for offsetof */

/* We use C99s 'bool' type, particularly for carries.  So we need to include
   stdbool.h if present, else define substitutes for bool, true and false */

#ifdef HAVE_STDBOOL_H
#  include <stdbool.h>
#else
#  ifndef HAVE__BOOL
#    define _Bool signed char;
#  endif
#  define bool _Bool;
#  define false 0;
#  define true 1;
#endif

/*
  A Deccoeff instance is stored internally in base LIMB_BASE =
  10**LIMB_DIGITS, as an array of limbs (least significant first).
  LIMB_DIGITS is 18 if 64-bit and 128-bit integer types are available, and 9
  if 32-bit and 64-bit integer types are available (as they should be on
  almost all modern machines); otherwise it's 4.

  (It's possible to squeeze 19 decimal digits into a 128-bit unsigned integer
  type, but that leaves little room for maneuver.  Using 18 decimal digits
  instead makes the basic multiplication algorithm significantly faster, by
  saving on the number of divisions necessary.)

  Here we define the various types and constants needed.

  limb_t should be an integer type that can hold any integer in the range [0,
  LIMB_BASE).  double_limb_t should be able to hold any integer in the range
  [0, LIMB_BASE*LIMB_BASE).  Type digit_limb_t is used for base conversion
  from base 2 to base 10 and back again and should be able to hold any
  integer in the range [0, LIMB_BASE * PyLong_BASE).

  BASEC_P/BASEC_Q is an upper bound for log(PyLong_BASE)/log(LIMB_BASE).
  BASECI_P/BASECI_Q is an upper bound for log(LIMB_BASE)/log(PyLong_BASE).
  These upper bounds are used for determining how much memory to allocate
  when doing base conversion.  BASEC_P * BASEC_Q and BASECI_P * BASECI_Q
  should fit comfortably in a Py_ssize_t.
*/

/* Rational approximations to log(10)/log(2), used for base conversion:
   485/146  = 3.321917808219...
   log2(10) = 3.321928094887...
   2136/643 = 3.321928460342...
*/

#define LOG2_10LP 485
#define LOG2_10LQ 146
#define LOG2_10UP 2136
#define LOG2_10UQ 643

#if defined(HAVE_STDINT_H)
#include <stdint.h>
#endif

#if defined(HAVE_INTTYPES_H)
#include <inttypes.h>
#endif

#if (defined(UINT64_MAX) || defined(uint64_t)) &&       \
    (defined(HAVE___UINT128_T))

/* if a 128-bit unsigned integer type is available, use a 64-bit limb with 18
   digits to a limb.  Should make this configurable, since it's quite possible
   that a 32-bit limb is faster, even on 64-bit machines. */

typedef uint64_t limb_t;
typedef __uint128_t double_limb_t;
typedef __uint128_t digit_limb_t;
#define LIMB_DIGITS 18
#define LIMB_MAX ((limb_t)999999999999999999) /* 10**LIMB_DIGITS - 1 */

#elif (defined(UINT32_MAX) || defined(uint32_t)) &&     \
    (defined(UINT64_MAX) || defined(uint64_t))

/* ... else use uint32_t for limb_t and uint64_t for double_limb_t if
   available, with 9 digits to a limb... */

typedef uint32_t limb_t;
typedef uint64_t double_limb_t;
typedef uint64_t digit_limb_t;
#define LIMB_DIGITS 9
#define LIMB_MAX ((limb_t)999999999)  /* 10**LIMB_DIGITS - 1 */

#else

/* ... else fall back to short and long, and make LIMB_DIGITS 4. */

typedef unsigned short limb_t;
typedef unsigned long double_limb_t;
typedef unsigned long digit_limb_t;
#define LIMB_DIGITS 4
#define LIMB_MAX ((limb_t)9999)

#endif

static limb_t powers_of_ten[LIMB_DIGITS];
#define LIMB_ZERO ((limb_t)0)
#define LIMB_ONE ((limb_t)1)

/*
  The rest of this file is organized in three parts.  First we have primitive
  operations on single limbs (arithmetic, shifts, etc.)---e.g. limb_adc,
  limb_fma, ....  This is followed by functions defining operations on arrays
  of limbs (limbs_add, limbs_multiply, ...); these functions contain the
  necessary mathematical algorithms, but know nothing about allocating memory,
  raising exceptions, converting types, etc.  Lastly we have functions
  providing the deccoeff operations.
*/


/*********************************
 * Primitive operations on limbs *
 *********************************/

/* The idea behind this section is that it should be easy to change the
   underlying representation of a limb without greatly affecting the rest of
   the code.  For example, on a 64-bit platform one might want to experiment
   with a 64-bit limb containing 18 or 19 digits per limb, or one might want
   to use alternative encodings like binary-coded decimal or densely packed
   decimal.  In that case, only this section should have to be changed; all
   following code uses only these primitive operations to operate on limbs.

   In particular, later code should not assume that limbs can be operated on
   using the usual arithmetic operators, or that they are comparable with < or
   == (some encodings may not be monotonic, or may have redundant encodings of
   the same integer, or may not even be encoded as a C integer type).
*/

#define LIMB_BASE (LIMB_MAX+(limb_t)1)

static void
limb_error(const char *msg)
{
    fprintf(stderr, "%s\n", msg);
    abort();
}

/* add with carry: compute a + b + c; put sum in *r and return new carry. */

static bool
limb_adc(limb_t *r, limb_t a, limb_t b, bool c)
{
    limb_t sum;
    sum = a + b + (c ? LIMB_ONE : LIMB_ZERO);
    if (sum >= LIMB_BASE) {
        *r = sum - LIMB_BASE;
        return true;
    }
    else {
        *r = sum;
        return false;
    }
}

/* subtract with borrow: compute a - (b + c); put result in *r and return the
   carry (true if a - (b+c) < 0, false otherwise). */

static bool
limb_sbb(limb_t *r, limb_t a, limb_t b, bool c)
{
    limb_t diff;
    diff = a - b - (c ? LIMB_ONE : LIMB_ZERO);
    if (diff > a) {
        *r = diff + LIMB_BASE;
        return true;
    }
    else {
        *r = diff;
        return false;
    }
}

/* multiply and add with two addends; this is useful for long multiplication.
   Store the low part of the result in *low, and return the high part. */

static limb_t
limb_fmaa(limb_t *low, limb_t a, limb_t b, limb_t c, limb_t d) {
    double_limb_t hilo;
    hilo = (double_limb_t)a * b + c + d;
    *low = (limb_t)(hilo%LIMB_BASE);
    return (limb_t)(hilo/LIMB_BASE);
}

/* division: divide high*BASE+low by c.  Return the quotient and stored
   remainder in *rem.  Requires that high < c. */

static limb_t
limb_div(limb_t *rem, limb_t high, limb_t low, limb_t c) {
    double_limb_t hilo;
    if (high >= c)
        limb_error("limb_div: invalid division");
    hilo = (double_limb_t)high * LIMB_BASE + low;
    *rem = (limb_t)(hilo%c);
    return (limb_t)(hilo/c);
}

/* determine whether the given limb is nonzero */

static bool
limb_bool(limb_t a)
{
    return a != 0;
}

/* The following two functions are used in base conversion.  Any value n in
   the range [0, LIMB_BASE * PyLong_BASE) can be written either in the form

   n = a * LIMB_BASE + b     'digit-limb'

   with 0 <= a < PyLong_BASE, 0 <= b < LIMB_BASE, or in the form

   n = c * PyLong_BASE + d    'limb-digit'

   with 0 <= c < LIMB_BASE, 0 <= d < PyLong_BASE.  digit_limb_swap converts
   from the first form to the second, limb_digit_swap does the reverse. */

static digit
digit_limb_swap(limb_t *c, digit a, limb_t b)
{
    digit_limb_t hilo;
    hilo = (digit_limb_t)a * LIMB_BASE + b;
    *c = (limb_t)(hilo >> PyLong_SHIFT);
    return (digit)(hilo & PyLong_MASK);
}

static limb_t
limb_digit_swap(digit *a, limb_t c, digit d)
{
    digit_limb_t hilo;
    hilo = ((digit_limb_t)c << PyLong_SHIFT) + d;
    *a = (digit)(hilo / LIMB_BASE);
    return (limb_t)(hilo % LIMB_BASE);
}

/* get a hash value from a limb */

static long
limb_hash(limb_t x) {
    return (long)x;
}

/* return index of most significant digit of given limb; result is undefined
   if the limb is 0 */

static Py_ssize_t
limb_msd(limb_t x) {
    Py_ssize_t i;
    if (x == 0)
        limb_error("limb_msd: zero argument");
    for (i=0; i < LIMB_DIGITS && powers_of_ten[i] <= x; i++);
    return i;
}

/* here's an attempt at faster multiplication; it's in this section because it
   depends on knowing the representation of a limb.  It saves on the number of
   divisions required by accumulating the sum of several partial products at a
   time.  The largest number of partials we can accumulate is 42 if
   LIMB_DIGITS == 4, 18 if LIMB_DIGITS == 9, and 340 if LIMB_DIGITS = 18. */

#if LIMB_DIGITS == 4
#define MAX_PARTIALS 42  /* floor(2^32 / 10^8) */
#elif LIMB_DIGITS == 9
#define MAX_PARTIALS 18  /* floor(2^64 / 10^18) */
#elif LIMB_DIGITS == 18
#define MAX_PARTIALS 340 /* floor(2^128 / 10^36) */
#else
#error "unrecognised value for LIMB_DIGITS"
#endif

/* res[0:a_size+b_size] := a*b, assuming b_size <= MIN(MAX_PARTIALS,
   a_size) */

static void
limbs_multiply_init(limb_t *res, const limb_t *a, Py_ssize_t a_size,
                const limb_t *b, Py_ssize_t b_size)
{
    double_limb_t acc = 0;
    Py_ssize_t j, k;
    assert(b_size <= MAX_PARTIALS && b_size <= a_size);
    for (k=0; k<b_size; k++) {
        for (j=0; j<=k; j++)
            acc += (double_limb_t)a[k-j]*b[j];
        res[k] = acc % LIMB_BASE;
        acc /= LIMB_BASE;
    }
    for (; k < a_size; k++) {
        for (j=0; j<b_size; j++)
            acc += (double_limb_t)a[k-j]*b[j];
        res[k] = acc % LIMB_BASE;
        acc /= LIMB_BASE;
    }
    for (k=0; k<b_size; k++) {
        for (j=k+1; j < b_size; j++)
            acc += (double_limb_t)a[k+a_size-j]*b[j];
        res[k+a_size] = acc % LIMB_BASE;
        acc /= LIMB_BASE;
    }
    assert(acc == 0);
}

/* res[0:a_size+b_size] := a*b + res[0:a_size], assuming b_size <=
   MIN(MAX_PARTIALS, a_size) */

static void
limbs_multiply_add(limb_t *res, const limb_t *a, Py_ssize_t a_size,
                const limb_t *b, Py_ssize_t b_size)
{
    double_limb_t acc = 0;
    Py_ssize_t j, k;
    assert(b_size <= MAX_PARTIALS && b_size <= a_size);
    for (k=0; k<b_size; k++) {
        acc += res[k];
        for (j=0; j<=k; j++)
            acc += (double_limb_t)a[k-j]*b[j];
        res[k] = acc % LIMB_BASE;
        acc /= LIMB_BASE;
    }
    for (; k < a_size; k++) {
        acc += res[k];
        for (j=0; j<b_size; j++)
            acc += (double_limb_t)a[k-j]*b[j];
        res[k] = acc % LIMB_BASE;
        acc /= LIMB_BASE;
    }
    for (k=0; k<b_size; k++) {
        for (j=k+1; j < b_size; j++)
            acc += (double_limb_t)a[k+a_size-j]*b[j];
        res[k+a_size] = acc % LIMB_BASE;
        acc /= LIMB_BASE;
    }
    assert(acc == 0);
}

/* res[0:a_size+b_size] := a * b */

static void
limbs_mul(limb_t *res, const limb_t *a, Py_ssize_t a_size,
                const limb_t *b, Py_ssize_t b_size)
{
    /* reduce to case where a_size >= b_size */
    if (a_size < b_size) {
        const limb_t *temp;
        Py_ssize_t temp_size;
        temp = a; a = b; b = temp;
        temp_size = a_size; a_size = b_size; b_size = temp_size;
    }

    assert(b_size <= a_size);
    if (b_size < MAX_PARTIALS)
        limbs_multiply_init(res, a, a_size, b, b_size);
    else {
        limbs_multiply_init(res, a, a_size, b, MAX_PARTIALS);
        b_size -= MAX_PARTIALS;
        b += MAX_PARTIALS;
        res += MAX_PARTIALS;
        while (b_size >= MAX_PARTIALS) {
            limbs_multiply_add(res, a, a_size, b, MAX_PARTIALS);
            b_size -= MAX_PARTIALS;
            b += MAX_PARTIALS;
            res += MAX_PARTIALS;
        }
        limbs_multiply_add(res, a, a_size, b, b_size);
    }
}

#undef LIMB_BASE


/*******************************
 * Derived operations on limbs *
 *******************************/

/* These limb operations are derived from the primitive operations above, and
   provided for convenience.  There is no need to change these if/when the
   representation of a limb_t changes. */

/* comparisons */

static int
limb_cmp(limb_t a, limb_t b)
{
    bool carry;
    limb_t diff;
    carry = limb_sbb(&diff, a, b, false);
    if (limb_bool(diff))
        return carry ? -1 : 1;
    else
        return 0;
}

/* a == b iff a - b is zero */

static bool
limb_eq(limb_t a, limb_t b)
{
    limb_t diff;
    limb_sbb(&diff, a, b, false);
    return !limb_bool(diff);
}

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

/* extract bottom n digits of a limb */

static limb_t
limb_mask(limb_t a, Py_ssize_t n)
{
    if (!(0 <= n && n <= LIMB_DIGITS))
        limb_error("limb_mask: invalid count");
    if (n < LIMB_DIGITS)
        limb_div(&a, LIMB_ZERO, a, powers_of_ten[n]);
    return a;
}

/* convert a character in the range '0' through '9' into a limb and back */

static limb_t
wdigit_to_limb(Py_UNICODE d)
{
    digit dummy;
    return limb_digit_swap(&dummy, LIMB_ZERO, (digit)(d - '0'));
}

static char
limb_to_digit(limb_t b)
{
    limb_t dummy;
    return '0' + (char)digit_limb_swap(&dummy, 0, b);
}

/* *res = a << n + b, b < 10**n.  Returns part shifted out. */

static limb_t
limb_lshift(limb_t *res, limb_t a, Py_ssize_t n, limb_t b) {
    if (!(0 <= n && n <= LIMB_DIGITS))
        limb_error("limb_lshift: invalid shift index");
    if (LIMB_DIGITS == n) {
        *res = b;
        return a;
    }
    else
        return limb_fmaa(res, a, powers_of_ten[n], b, LIMB_ZERO);
}

/* *res = (a + b*LIMB_BASE) >> n, b < 10**n.  Returns part shifted out. */

static limb_t
limb_rshift(limb_t *res, limb_t a, Py_ssize_t n, limb_t b) {
    limb_t rem;
    if (!(0 <= n && n <= LIMB_DIGITS))
        limb_error("limb_lshift: invalid shift index");
    if (LIMB_DIGITS == n) {
        *res = b;
        return a;
    }
    else {
        *res = limb_div(&rem, b, a, powers_of_ten[n]);
        return rem;
    }
}

/* retrieve the value of a particular digit, as a limb_t */

static limb_t
limb_getdigit(limb_t x, Py_ssize_t n)
{
    limb_t q, dummy;
    if (!(0 <= n && n < LIMB_DIGITS))
        limb_error("invalid digit in limb_getdigit");
    q = limb_div(&dummy, LIMB_ZERO, x, powers_of_ten[n]);
    return limb_mask(q, 1);
}

/*********************************
 * Arithmetic on arrays of limbs *
 *********************************/

/* Low-level functions for operating on arrays of limbs.  These functions
   don't take care of memory allocation; they assume that sufficient space is
   provided for their results. */

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

/* increment a_size-limb number a if carry is true, else just copy it; gives
   a_size-limb result and returns a carry */

static bool
limbs_incc(limb_t *res, const limb_t *a, Py_ssize_t a_size, bool carry)
{
    Py_ssize_t i;
    for (i=0; i < a_size; i++)
        carry = limb_adc(res+i, a[i], LIMB_ZERO, carry);
    return carry;
}

/* decrement a_size-limb number a if carry is true, else just copy it; gives
   a_size-limb result and returns a carry */

static bool
limbs_decc(limb_t *res, const limb_t *a, Py_ssize_t a_size, bool carry)
{
    Py_ssize_t i;
    for (i=0; i < a_size; i++)
        carry = limb_sbb(res+i, a[i], LIMB_ZERO, carry);
    return carry;
}

/* subtract limbs from zero, with borrow.  Equivalent to the nines' complement
   if carry == true on input, and the ten's complement otherwise. */

static bool
limbs_cmpl(limb_t *res, const limb_t *a, Py_ssize_t a_size, bool carry)
{
    Py_ssize_t i;
    for (i=0; i < a_size; i++)
        carry = limb_sbb(res+i, LIMB_ZERO, a[i], carry);
    return carry;
}

/* copy limbs */

static void
limbs_copy(limb_t *res, const limb_t *a, Py_ssize_t a_size)
{
    Py_ssize_t i;
    for (i=0; i < a_size; i++)
        res[i] = a[i];
}

/* write zero to a set of limbs */

static void
limbs_zero(limb_t *res, Py_ssize_t n)
{
    Py_ssize_t i;
    for (i=0; i < n; i++)
        res[i] = LIMB_ZERO;
}

/* take the absolute value of the difference between n-limb numbers a and b,
   returning true if a < b, false if a >= b. */

static bool
limbs_diff(limb_t *res, const limb_t *a, const limb_t *b, Py_ssize_t n)
{
    Py_ssize_t i;
    bool carry, neg;

    /* reduce to case where top limbs differ */
    while (n > 0 && limb_eq(a[n-1], b[n-1])) {
        res[n-1] = LIMB_ZERO;
        n--;
    }
    if (n == 0)
        return false;

    /* now reduce to case a > b */
    neg = limb_lt(a[n-1], b[n-1]);
    if (neg) {
        const limb_t *temp;
        temp = a; a = b; b = temp;
    }

    /* subtract */
    carry = false;
    for (i=0; i < n; i++)
        carry = limb_sbb(res+i, a[i], b[i], carry);
    assert(!carry);
    return neg;
}

/* multiply a_size-limb number a by single limb x, getting a_size-limb result
   res and returning the high limb. */

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

/* divide a_size-limb number a by single limb x, giving a_size-limb quotient
   res and returning the (single limb) remainder */

static limb_t
limbs_div1(limb_t *res, const limb_t *a, Py_ssize_t a_size,
           limb_t high, limb_t x)
{
    Py_ssize_t i;
    for (i = a_size-1; i >= 0; i--)
        res[i] = limb_div(&high, high, a[i], x);
    return high;
}

/* divide a by b, giving an (a_size-b_size+1)-limb quotient and b_size-limb
   remainder.  Assumes that the top limb of b is nonzero and that a_size >=
   b_size.  w provides a_size+b_size+1 limbs of workspace. */

static void
limbs_div(limb_t *quot, limb_t *rem, const limb_t *a, Py_ssize_t a_size,
          const limb_t *b, Py_ssize_t b_size, limb_t *w)
{
    limb_t scale, top, a_top, b_top, q, dummy;
    limb_t *aa, *bb;
    bool carry;
    Py_ssize_t j;

    /* top limb of b should be nonzero; a should contain at least as many
       limbs as b */
    assert(a_size >= b_size && b_size > 0 && limb_bool(b[b_size-1]));

    /* compute scale factor for normalization: floor(LIMB_BASE /
       (b_top+1)) */
    carry = limb_adc(&scale, b[b_size-1], LIMB_ONE, false);
    if (carry)
        scale = LIMB_ONE;
    else
        scale = limb_div(&dummy, LIMB_ONE, LIMB_ZERO, scale);

    /* scale a and b */
    top = limbs_mul1(w, b, b_size, scale);
    bb = w;
    assert(!limb_bool(top));

    top = limbs_mul1(w+b_size, a, a_size, scale);
    aa = w+b_size;

    /* catch most cases where quotient only needs a_size-b_size limbs */
    if (!limb_bool(top) && limb_lt(aa[a_size-1], bb[b_size-1]))
        quot[a_size-b_size] = LIMB_ZERO;
    else {
        aa[a_size] = top;
        a_size++;
    }

    b_top = bb[b_size-1];
    aa += a_size-b_size;
    for (j = a_size-b_size-1; j >= 0; j--) {
        aa--;
        a_top = aa[b_size];
        assert(limb_le(a_top, b_top));
        /* quotient q = aa / bb; may be overestimate */
        if (limb_lt(a_top, b_top))
            q = limb_div(&dummy, a_top, aa[b_size-1], b_top);
        else
            q = LIMB_MAX;
        /* compute bottom b_size limbs of aa[j:] - q*bb */
        top = limbs_mul1(rem, bb, b_size, q);
        carry = limbs_sub(aa, aa, rem, b_size);
        carry = limb_adc(&top, top, LIMB_ZERO, carry);
        assert(!carry);
        assert(limb_le(a_top, top));
        /* correct if necessary */
        while (limb_lt(a_top, top)) {
            carry = limbs_add(aa, aa, bb, b_size);
            carry = limb_adc(&a_top, a_top, LIMB_ZERO, carry);
            assert(!carry);
            carry = limb_sbb(&q, q, LIMB_ONE, false);
            assert(!carry);
        }
        quot[j] = q;
    }
    top = limbs_div1(rem, aa, b_size, LIMB_ZERO, scale);
    assert(!limb_bool(top));
}

/* shift a_size-limb number left n digits (shifting zeros in); i.e., multiply
   by 10**n.  Fills the first a_size + ceiling(n/LIMB_DIGITS) limbs of res. */

static void
limbs_lshift(limb_t *res, const limb_t *a, Py_ssize_t a_size, Py_ssize_t n)
{
    Py_ssize_t n_limbs, n_digits, i;
    limb_t high;
    assert(n >= 0);
    n_limbs = n / LIMB_DIGITS;
    n_digits = n % LIMB_DIGITS;
    for (i = 0; i < n_limbs; i++)
        res[i] = LIMB_ZERO;
    high = limbs_mul1(res+n_limbs, a, a_size, powers_of_ten[n_digits]);
    assert(limb_lt(high, powers_of_ten[n_digits]));
    if (n_digits != 0)
        res[n_limbs + a_size] = high;
}

/* shift a_size-limb number a right n digits (not limbs!).  i.e., divide by
   10**n.  n should satisfy 0 <= n <= LIMB_DIGITS*a_size, and the result res
   has a_size - n/LIMB_DIGITS limbs. */

static void
limbs_rshift(limb_t *res, const limb_t *a, Py_ssize_t a_size, Py_ssize_t n)
{
    Py_ssize_t n_limbs, n_digits;
    assert(0 <= n && n <= a_size*LIMB_DIGITS);
    n_limbs = n / LIMB_DIGITS;
    n_digits = n % LIMB_DIGITS;
    limbs_div1(res, a+n_limbs, a_size-n_limbs, LIMB_ZERO,
               powers_of_ten[n_digits]);
}

/* get slice a[m:n] of (the decimal digits of) an integer a, from digit m up
   to (but not including) digit n, giving a result res with
   ceiling((n-m)/LIMB_DIGITS) limbs.  Assumes that 0 <= m < n and that a has
   at least 1+(n-1)/LIMB_DIGITS limbs. */

static void
limbs_slice(limb_t *res, const limb_t *a, Py_ssize_t m, Py_ssize_t n)
{
    Py_ssize_t mlimbs, mdigits, reslimbs, resdigits, diff;
    limb_t high;
    mlimbs = m / LIMB_DIGITS;
    mdigits = m % LIMB_DIGITS;
    reslimbs = (n-1-m) / LIMB_DIGITS;
    resdigits = (n-1-m) % LIMB_DIGITS;
    diff = (mdigits + resdigits - LIMB_DIGITS);
    if (diff < 0) {
        limbs_div1(res, a + mlimbs, reslimbs + 1,
                   LIMB_ZERO, powers_of_ten[mdigits]);
        res[reslimbs] = limb_mask(res[reslimbs], resdigits+1);
    }
    else {
        high = limb_mask(a[mlimbs + reslimbs + 1], diff+1);
        limbs_div1(res, a + mlimbs, reslimbs + 1,
                   high, powers_of_ten[mdigits]);
    }
}

/* get a particular digit from an array of limbs.  assumes that 0 <= n <
   a_size * LIMB_DIGITS */

static limb_t
limbs_getdigit(limb_t *a, Py_ssize_t n)
{
    return limb_getdigit(a[n/LIMB_DIGITS], n%LIMB_DIGITS);
}

/* Conversion to and from strings */

/* Fill a character array from an array of limbs */

static void
limbs_as_unicode(Py_UNICODE *s, Py_ssize_t s_len, const limb_t *a)
{
    Py_UNICODE *s_store;
    limb_t limb;
    Py_ssize_t nlimbs, ndigits, i, j;
    if (s_len == 0)
        return;

    /* s_len == nlimbs*LIMB_DIGITS + ndigits, 0 < ndigits <= LIMB_DIGITS */
    nlimbs = (s_len-1)/LIMB_DIGITS;
    ndigits = (s_len-1)%LIMB_DIGITS + 1;

    /* fill in digits from right to left */
    s_store = s;
    s += s_len;
    for (j=0; j < nlimbs; j++) {
        limb = a[j];
        for (i=0; i < LIMB_DIGITS; i++)
            *--s = limb_to_digit(limb_rshift(&limb, limb, 1, LIMB_ZERO));
    }
    /* most significant limb */
    limb = a[nlimbs];
    for (i=0; i < ndigits; i++)
        *--s = limb_to_digit(limb_rshift(&limb, limb, 1, LIMB_ZERO));

    assert(s == s_store);
}


/* Convert a character array to an array of limbs.  Returns false on success,
   true if any of the characters was not a valid digit; in the latter case,
   the contents of a are undefined. a should provide ceiling(slen /
   LIMB_DIGITS) limbs. */

static bool
limbs_from_unicode(limb_t *a, const Py_UNICODE *s, Py_ssize_t s_len)
{
    Py_ssize_t i, k, digits_in_limb;
    limb_t acc;
    Py_UNICODE c;

    k = (s_len+LIMB_DIGITS-1) / LIMB_DIGITS;
    digits_in_limb = (s_len+LIMB_DIGITS-1) % LIMB_DIGITS + 1;
    acc = LIMB_ZERO;
    for (i = 0; i < s_len; i++) {
        c = s[i];
        if (c < '0' || c > '9')
            return true;
        limb_lshift(&acc, acc, 1, wdigit_to_limb(c));
        digits_in_limb--;
        if (digits_in_limb == 0) {
            digits_in_limb = LIMB_DIGITS;
            a[--k] = acc;
            acc = LIMB_ZERO;
        }
    }
    assert(digits_in_limb == LIMB_DIGITS);
    return false;
}

/*
   variant of limbs_from_string that accepts a string where all but one of the
   characters is a decimal digit, and ignores the specified character.
   (Typically, the character to be ignored is a decimal point.)

   On entry, s is a pointer to a character array of length s_len + 1, 0 <=
   int_len <= s_len, s[0] through s[int_len-1] and s[int_len+1] through
   s[s_len] contain decimal digits.  The result is then identical to
   limbs_from_string applied to the concatenation of s[:int_len] and
   s[int_len+1:].  No checking is performed, and there is no return value.

   This function is used when parsing _Decimal instances.
 */

static void
limbs_from_pointed_unicode(limb_t *a, const Py_UNICODE *s, Py_ssize_t s_len,
                          Py_ssize_t int_len)
{
    Py_ssize_t i, k, digits_in_limb;
    limb_t acc;
    Py_UNICODE c;

#define GET_DIGIT(i) ((i) < int_len ? s[(i)] : s[(i)+1])

    k = (s_len+LIMB_DIGITS-1) / LIMB_DIGITS;
    digits_in_limb = (s_len+LIMB_DIGITS-1) % LIMB_DIGITS + 1;
    acc = LIMB_ZERO;
    for (i = 0; i < s_len; i++) {
        c = GET_DIGIT(i);
        limb_lshift(&acc, acc, 1, wdigit_to_limb(c));
        digits_in_limb--;
        if (digits_in_limb == 0) {
            digits_in_limb = LIMB_DIGITS;
            a[--k] = acc;
            acc = LIMB_ZERO;
        }
    }
    assert(digits_in_limb == LIMB_DIGITS);
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
            a[i] = limb_digit_swap(&high, a[i], high);
        while (high != 0)
            a[a_size++] = limb_digit_swap(&high, LIMB_ZERO, high);
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
            b[j] = digit_limb_swap(&high, b[j], high);
        while (limb_bool(high))
            b[b_size++] = digit_limb_swap(&high, 0, high);
    }
    return b_size;
}

/* print a limb array; debugging aid */

static void
limbs_printf(const char *name, const limb_t *a, Py_ssize_t a_size)
{
    Py_ssize_t i;
    printf("%s = [", name);

    for (i=0; i<a_size; i++) {
        if (i != 0)
            printf(", ");
        printf("%09u", a[i]);
    }
    printf("]\n");
}

/****************************
 * Karatsuba multiplication *
 ****************************/

/* Note that for a balanced multiplication, where you don't want to
   keep the original multiplicands, one can do the multiplication
   without any additional workspace necessary:

   to multiply a0 a1 * b0 b1, with result in z0 z1 z2 z3
   (1) a0 - a1 -> z2
   (2) b0 - b1 -> z3
   (3) a0 * b0 -> z0 z1
   (4) b1 -> a0
   (5) z2 * z3 -> b0 b1
   (6) a0 * a1 -> z2 z3
   (7)
*/


/* Unbalanced karatsuba multiplication: multiply a and b, putting result (of
   size a_size + b_size) in res.  w provides workspace.

   Preconditions: 1 <= k < a_size <= b_size <= 2*k, and k is a
   power of 2.

   Algorithm
   ---------
   Write m for a_size, n for b_size, B for LIMB_BASE.  Then we're assuming
   that

   k < m <= 2*k, k < n <= 2*k

   and so in particular, 0 < m-k <= k, 0 < n-k <= k.

   Let a0, a1 = a[:k], a[k:]
   Let b0, b1 = b[:n-k], b[n-k:]
   Then a = a0 + a1*B**k, b = b0 + b1*B**(n-k)

   Let c = a0 - a1, d = b1 - b0*B**(2k-n)
   Then

   ab = a0*b0 + (c*d + a1*b1 + a0*b0*B**(2k-n))*B**(n-k) + a1*b1*B**n.

   +---+---+---+
   +   a0  |a1 +
   +---+---+---+
   0  m-k  k   m

   +-----+-+-----+
   | b0  | b1    |
   +-----+-+-----+
   0   n-k k     n

   Sizes:
   a0    k
   a1    m-k
   b0    n-k
   b1    k
   c     k
   d     k
   a0*b0 n
   a1*b1 m
   c*d   2*k
*/


/* notes:

   (1) there's no real need to pass w_size around, but it helps guard
   against errors
   (2) slice notation in the comments below is as in Python:
   res[2n:4n] means the number represented by res[2*n] through
   res[4*n-1], inclusive.
*/

/* decide how to handle a particular size of multiplication, and
   dispatch to the appropriate function to do the computation */

/* strategy: wlog suppose 1 <= m <= n.

   (1) if m == 1, do a basecase multiplication.

   Otherwise, let k be the largest power of 2 strictly less than m
   (so k < m <= 2*k).  Then:

   (2) If n <= 2*k, do a Karatsuba multiplication.  Else

   (3) Otherwise, split n as n[:2k] and n[2k:], and do an m*(2k) Karatsuba
   multiplication and a m * (n-2k) multiplication (via a recursive call to
   limbs_mul_dispatch).  Add the results.

   Question: how much workspace do we need for this?

   Write W(m, n) for the amount of workspace needed; assume 1 <= m <= n.

   Case 0: W(1, n).  No workspace required.

   Case 1: W(k, k), k a power of 2:
   W(k, k) = 2k+max(1, W(k/2, k/2)); W(1, 1) = 0.  So W(k, k) = 2k-1.

   Case 2: W(m, n),  k < m <= n <= 2k, k a power of 2:
   W(m, n) = 4k-2, assuming that W(m, k) <= 2k-2 for all m < k...

   Case 3: W(m, k),  m <= k, k a power of 2:
   Case 3a: if k/2 < m <= k then by Case 2, W(m, k) = 2k-2.
   Case 3b: if k/4 < m <= k/2 we'll end up doing 2 m-by-k/2 Karatsubas,
   needing a total of k-2 (workspace for Karatsubas)
   plus (m+k/2) for storage of result of 2nd Karatsuba
   giving 3k/2 -2 + m <= 2k-2.
   ... and so on.  In all cases, 2k-2 should be sufficient.

   Case 4: W(m, n), k < m <= 2k < n, k a power of 2.  Then
   we do a number of m-by-2k multiplications;  for the second
   and subsequent of these, we need to store the result in
   the workspace as well, so workspace required is at most
   W(m, 2k) + m + 2k = 6k-2 + m <= 8k-2.

   Conclusion: let k < m <= 2k, l < n <= 2l, k and l powers of 2.
   Then for l == k, 4k-2 limbs are enough
   for l >= 2k, 8k-2 limbs are enough
   In all cases, 4l-2 limbs are enough.
*/

#define KARATSUBA_CUTOFF 72

/* limbs_kmul and limbs_mul_dispatch call each other recursively,
   so we need a forward declaration */

static void
limbs_mul_dispatch(limb_t *res, const limb_t *a, Py_ssize_t a_size,
                   const limb_t *b, Py_ssize_t b_size,
                   limb_t *w, Py_ssize_t w_size);

/* Karatsuba multiplication.  On input, k is an integer satisfying k < a_size
   <= 2*k and k < b_size <= 2*k.  w provides workspace of size w_size.  The
   amount of workspace required is 2k + max(1, W), where W is the workspace
   required by the recursive calls.*/

static void
limbs_kmul(limb_t *res, const limb_t *a, Py_ssize_t a_size,
           const limb_t *b, Py_ssize_t b_size, Py_ssize_t k,
           limb_t *w, Py_ssize_t w_size)
{
    bool carry, sign;

    /* check preconditions */
    assert(k < a_size && a_size <= 2*k && k < b_size && b_size <= 2*k);

    /* abs(a0 - a1) */
    limbs_copy(res, a+k, a_size-k);
    limbs_zero(res+a_size-k, 2*k-a_size);
    sign = limbs_diff(res, a, res, k);

    /* abs(b1 - b0*B**(2k-n)) */
    limbs_zero(res+k, 2*k-b_size);
    limbs_copy(res+k+(2*k-b_size), b, b_size-k);
    sign ^= limbs_diff(res+k, b+b_size-k, res+k, k);

    /* need 2k+1 limbs to store central product */
    assert(w_size >= 2*k+1);
    limbs_mul_dispatch(w, res, k, res+k, k, w+2*k, w_size-2*k);
    /* products a0*b0 and a1*b1 are placed directly into res */
    limbs_mul_dispatch(res, b, b_size-k, a, k, w+2*k, w_size-2*k);
    limbs_mul_dispatch(res+b_size, a+k, a_size-k, b+b_size-k, k,
                       w+2*k, w_size-2*k);
    if (sign) {
        carry = limbs_sub(w, res+b_size, w, a_size);
        carry = limbs_cmpl(w+a_size, w+a_size, 2*k-a_size, carry);
    }
    else {
        carry = limbs_add(w, res+b_size, w, a_size);
        carry = limbs_incc(w+a_size, w+a_size, 2*k-a_size, carry);
    }
    carry ^= limbs_add(w+2*k-b_size, w+2*k-b_size, res, b_size);
    w[2*k] = carry ? LIMB_ONE : LIMB_ZERO;
    /* add central term in */
    carry = limbs_add(res+b_size-k, res+b_size-k, w, 2*k+1);
    carry = limbs_incc(res+b_size+k+1, res+b_size+k+1, a_size-k-1, carry);
    assert(!carry);
}

/* decide how to handle a particular size of multiplication, and
   dispatch to the appropriate function to do the computation.
   Assumes that 1 <= a_size <= b_size on input. */

static void
limbs_mul_dispatch(limb_t *res, const limb_t *a, Py_ssize_t a_size,
                   const limb_t *b, Py_ssize_t b_size,
                   limb_t *w, Py_ssize_t w_size)
{
    Py_ssize_t k;
    bool carry;
    assert(1 <= a_size && a_size <= b_size);

    if (a_size <= KARATSUBA_CUTOFF) {
        /* basecase multiplication */
        limbs_mul(res, a, a_size, b, b_size);
        return;
    }

    /* find k, the smallest number of the form 2**i * KARATSUBA_CUTOFF
       such that a_size <= k. */
    k = KARATSUBA_CUTOFF;
    while (k < a_size)
        k *= 2;
    assert(k/2 < a_size && a_size <= k);
    /* if sizes are balanced, do a single Karatsuba multiplication */
    if (b_size <= k) {
        limbs_kmul(res, a, a_size, b, b_size, k/2, w, w_size);
        return;
    }
    /* otherwise, while b_size > k, split off a k-by-a_size chunk and
       compute using Karatsuba.  The result from the first chunk goes
       straight into res; results from subsequent chunks go into the
       workspace and are then added into res */
    limbs_kmul(res, a, a_size, b, k, k/2, w, w_size);
    b += k; b_size -= k; res += k;
    while (b_size > k) {
        assert(w_size >= a_size+k);
        limbs_kmul(w, a, a_size, b, k, k/2,
                   w+a_size+k, w_size-a_size-k);
        carry = limbs_add(res, res, w, a_size);
        carry = limbs_incc(res+a_size, w+a_size, k, carry);
        assert(!carry);
        b += k; b_size -= k; res += k;
    }
    /* last chunk */
    assert(w_size >= a_size+b_size);
    if (b_size > k/2)
        limbs_kmul(w, a, a_size, b, b_size, k/2,
                   w+a_size+b_size, w_size-a_size-b_size);
    else
        limbs_mul_dispatch(w, b, b_size, a, a_size,
                           w+a_size+b_size, w_size-a_size-b_size);
    carry = limbs_add(res, res, w, a_size);
    carry = limbs_incc(res+a_size, w+a_size, b_size, carry);
    assert(!carry);
}

/**************************
 * deccoeff : definitions *
 **************************/

#define MODULE_NAME "deccoeff"
#define CLASS_NAME "Deccoeff"

/* We place an upper bound MAX_DIGITS on the number of decimal digits (*not*
   the number of limbs) in a deccoeff.  MAX_DIGITS should fit into a
   Py_ssize_t, so that length and indexing always make sense.  Here we take
   MAX_DIGITS = 10**9; then assuming that a Py_ssize_t is at least 32 bits
   we're safe from overflow even when adding two digit counts. */

#define MAX_DIGITS 1000000000

/* We really want ob_limbs[0] in the definition of deccoeff, but that's not
   standards-compliant C.  Use offsetof to get tp_basicsize. */

typedef struct {
    PyObject_VAR_HEAD
    limb_t ob_limbs[1];
} deccoeff;

#define DECCOEFF_ITEMSIZE sizeof(limb_t)
#define DECCOEFF_BASICSIZE (offsetof(deccoeff, ob_limbs))

static PyTypeObject deccoeff_DeccoeffType;

/* allocate a new decimal integer with 'size' limbs */

static deccoeff *
_deccoeff_new(Py_ssize_t size)
{
    /* XXX check for overflow */
    return PyObject_NEW_VAR(deccoeff, &deccoeff_DeccoeffType, size);
}

/* normalize a deccoeff, by setting the size correctly. */

static deccoeff *
deccoeff_normalize(deccoeff *v)
{
    Py_ssize_t v_size;
    v_size = Py_SIZE(v);
    while (v_size > 0 && !limb_bool(v->ob_limbs[v_size-1]))
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
        topdigits = limb_msd(v->ob_limbs[v_size-1]);
        small = topdigits <= (MAX_DIGITS-1)%LIMB_DIGITS+1;
    }
    else
        small = false;

    if (small)
        return v;
    Py_DECREF(v);
    PyErr_SetString(PyExc_OverflowError,
                    CLASS_NAME " instance has too many digits");
    return NULL;
}

static deccoeff *
_deccoeff_copy(deccoeff *a)
{
    Py_ssize_t a_size, i;
    deccoeff *z;

    a_size = Py_SIZE(a);
    z = _deccoeff_new(a_size);
    if (z != NULL)
        for (i=0; i < a_size; i++)
            z->ob_limbs[i] = a->ob_limbs[i];
    return z;
}

/* return zero */

static deccoeff *
deccoeff_zero(void)
{
    return _deccoeff_new(0);
}

/* return one */

static deccoeff *
deccoeff_one(void)
{
    deccoeff *z;
    z = _deccoeff_new(1);
    if (z != NULL)
        z->ob_limbs[0] = LIMB_ONE;
    return z;
}

static deccoeff *
_deccoeff_from_unicode_and_size(const Py_UNICODE *s, Py_ssize_t s_len) {
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

    invalid = limbs_from_unicode(z->ob_limbs, s, s_len);
    if (invalid) {
        Py_DECREF(z);
        PyErr_SetString(PyExc_ValueError,
                        "nondigit character in input");
        return NULL;
    }
    return deccoeff_normalize(z);
}

static deccoeff *
_deccoeff_from_pointed_unicode_and_size(const Py_UNICODE *s, Py_ssize_t s_len,
                                      Py_ssize_t int_len) {
    Py_ssize_t z_size;
    deccoeff *z;

    if (s_len > MAX_DIGITS) {
        PyErr_SetString(PyExc_OverflowError,
                        "too many digits");
        return NULL;
    }

    z_size = (s_len + LIMB_DIGITS - 1) / LIMB_DIGITS;
    z = _deccoeff_new(z_size);
    if (z == NULL)
        return NULL;

    limbs_from_pointed_unicode(z->ob_limbs, s, s_len, int_len);
    return deccoeff_normalize(z);
}



/***************************
 * Arithmetic on deccoeffs *
 ***************************/

/* General rules: if the result of any arithmetic operation falls
   outside the range [0, 10**MAX_DIGITS) then OverflowError is raised.
   Results are always normalized. */

/* determine whether another Python object is compatible with deccoeff, in the
   sense that it can be used in mixed-type arithmetic with deccoeff */

static bool
compatible_with_deccoeff(PyObject *v)
{
    return (v->ob_type == &deccoeff_DeccoeffType || PyIndex_Check(v));
}

/* convert an arbitrary PyObject to a deccoeff.  If conversion is implemented
   for this type, return false and put the result of the attempted conversion
   (which may be NULL on failure) in *z.   Otherwise, return true. */

static deccoeff *deccoeff_from_PyLong(PyLongObject *a);

static deccoeff *
convert_to_deccoeff(PyObject *v)
{
    PyLongObject *w;
    deccoeff *z;

    if (v->ob_type == &deccoeff_DeccoeffType) {
        Py_INCREF(v);
        return (deccoeff *)v;
    }
    else if (PyIndex_Check(v)) {
        w = (PyLongObject *)PyNumber_Index(v);
        if (w == NULL)
            return NULL;
        z = deccoeff_from_PyLong(w);
        Py_DECREF(w);
        return z;
    }
    else {
        PyErr_SetString(PyExc_TypeError,
                        "invalid type in deccoeff coercion");
        return NULL;
    }
}

#define DECCOEFF_WRAP_BINOP(PO_func, DC_func)                        \
    static PyObject *                                                \
    PO_func(PyObject *v, PyObject *w)                                \
    {                                                                \
        deccoeff *a, *b;                                        \
        PyObject *z = NULL;                                        \
        if (!compatible_with_deccoeff(v)) {                        \
            Py_INCREF(Py_NotImplemented);                        \
            z = Py_NotImplemented;                                \
        }                                                        \
        else if ((a = convert_to_deccoeff(v)) != NULL) {        \
            if (!compatible_with_deccoeff(w)) {                 \
                Py_INCREF(Py_NotImplemented);                   \
                z = Py_NotImplemented;                          \
            }                                                   \
            else if ((b = convert_to_deccoeff(w)) != NULL) {        \
                z = (PyObject *)(DC_func(a, b));                \
                Py_DECREF(b);                                   \
            }                                                   \
            Py_DECREF(a);                                        \
        }                                                        \
        return z;                                                \
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
    Py_ssize_t a_size, b_size, w_size;
    deccoeff *z, *w;

    a_size = Py_SIZE(a);
    b_size = Py_SIZE(b);
    if (a_size == 0 || b_size == 0)
        return deccoeff_zero();

    z = _deccoeff_new(a_size + b_size);
    if (z == NULL)
        return NULL;

    /* try out Karatsuba multiplication */

    /* swap so that a_size <= b_size */
    if (a_size > b_size) {
        Py_ssize_t temp_size;
        deccoeff *temp;
        temp = a; a = b; b = temp;
        temp_size = a_size; a_size = b_size; b_size = temp_size;
    }

    /* figure out how much workspace we need */
    /* find next power of 2 after a_size */
    w_size = 1;
    while (w_size < b_size)
        w_size *= 2;

    w_size = 2*w_size - 1;
    w = _deccoeff_new(w_size);
    if (w == NULL) {
        Py_DECREF(z);
        return NULL;
    }
    limbs_mul_dispatch(z->ob_limbs,
                       a->ob_limbs, a_size,
                       b->ob_limbs, b_size,
                       w->ob_limbs, w_size);

    Py_DECREF(w);
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
                                      a_size, LIMB_ZERO, b->ob_limbs[0]);
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

/* a*b % c;  assumes that a < c and b < c. */

static deccoeff *
_deccoeff_multiply_and_reduce(deccoeff *a, deccoeff *b, deccoeff *c)
{
    Py_ssize_t a_size, b_size;
    deccoeff *z, *w;

    a_size = Py_SIZE(a);
    b_size = Py_SIZE(b);
    z = _deccoeff_new(a_size + b_size);
    if (z == NULL)
        return NULL;
    limbs_mul(z->ob_limbs, a->ob_limbs, a_size, b->ob_limbs, b_size);
    /* w = z % c */
    w = _deccoeff_remainder(z, c);
    Py_DECREF(z);
    return w;
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


/* this version of power is naive and slow */

static deccoeff *
_deccoeff_power(deccoeff *a, deccoeff *bb, deccoeff *c)
{
    deccoeff *acc=NULL, *temp, *apow, *b;
    limb_t lowbit, *b_limbs;
    Py_ssize_t b_size;
    /* make a copy b of bb, which we'll modify in-place */
    b = _deccoeff_copy(bb);
    if (b == NULL)
        goto fail0;
    b_size = Py_SIZE(b);
    b_limbs = b->ob_limbs;
    /* apow keeps a power of a; starts as a (or a % c) */
    if (c == NULL) {
        apow = a;
        Py_INCREF(apow);
    }
    else
        apow = _deccoeff_remainder(a, c);
    if (apow == NULL)
        goto fail1;
    /* result accumulates in acc, which starts out as 1 (or 0 if c == 1) */
    acc = deccoeff_one();
    if (acc != NULL && c != NULL) {
        /* acc %= c */
        temp = _deccoeff_remainder(acc, c);
        Py_DECREF(acc);
        acc = temp;
    }
    if (acc == NULL || b_size == 0)
        goto fail2;
    while (true) {
        /* invariant quantity: apow**b*acc == a**bb. */
        lowbit = limbs_div1(b_limbs, b_limbs, b_size, LIMB_ZERO, 2);
        if (!limb_bool(b_limbs[b_size-1]))
            b_size--;
        if (limb_bool(lowbit)) {
            /* acc *= apow */
            if (c == NULL)
                temp = _deccoeff_multiply(apow, acc);
            else
                temp = _deccoeff_multiply_and_reduce(apow,
                                                     acc, c);
            Py_DECREF(acc);
            acc = temp;
        }
        if (acc == NULL || b_size == 0)
            break;

        /* apow *= apow */
        if (c == NULL)
            temp = _deccoeff_multiply(apow, apow);
        else
            temp = _deccoeff_multiply_and_reduce(apow, apow, c);
        Py_DECREF(apow);
        apow = temp;
        if (apow == NULL) {
            Py_DECREF(acc);
            acc = NULL;
            goto fail1;
        }
    }
  fail2:
    Py_DECREF(apow);
  fail1:
    Py_DECREF(b);
  fail0:
    return acc;
}

/* wrap _deccoeff_power */

static PyObject *
deccoeff_power(PyObject *v, PyObject *w, PyObject *x)
{
    deccoeff *a, *b, *c;
    PyObject *z;
    if (!compatible_with_deccoeff(v)) {
        Py_INCREF(Py_NotImplemented);
        return Py_NotImplemented;
    }
    a = convert_to_deccoeff(v);
    if (a == NULL)
        return NULL;

    if (!compatible_with_deccoeff(w)) {
        Py_DECREF(a);
        Py_INCREF(Py_NotImplemented);
        return Py_NotImplemented;
    }
    b = convert_to_deccoeff(w);
    if (b == NULL) {
        Py_DECREF(a);
        return NULL;
    }

    if (x == Py_None)
        c = NULL;
    else if (compatible_with_deccoeff(x)) {
        c = convert_to_deccoeff(x);
        if (c == NULL) {
            Py_DECREF(b);
            Py_DECREF(a);
            return NULL;
        }
    }
    else {
        Py_DECREF(b);
        Py_DECREF(a);
        Py_INCREF(Py_NotImplemented);
        return Py_NotImplemented;
    }

    z = (PyObject *)_deccoeff_power(a, b, c);

    Py_XDECREF(c);
    Py_DECREF(b);
    Py_DECREF(a);
    return z;
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
                        CLASS_NAME " instance has too many digits");
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
    else {
        start = PyNumber_AsSsize_t(b, NULL);
        if (start == -1 && PyErr_Occurred())
            return NULL;
        else if (start < 0) {
            PyErr_SetString(PyExc_ValueError,
                            "index should be nonnegative");
            return NULL;
        }
        if (start >= defstop)
            return deccoeff_zero();

        z = _deccoeff_new(1);
        if (z == NULL)
            return NULL;
        z->ob_limbs[0] = limbs_getdigit(a->ob_limbs, start);
        return deccoeff_normalize(z);
    }
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
    int c;
    Py_ssize_t a_size, b_size, i;
    a_size = Py_SIZE(a);
    b_size = Py_SIZE(b);
    if (a_size != b_size)
        return a_size < b_size ? -1 : 1;
    for (i = a_size-1; i >= 0; i--) {
        c = limb_cmp(a->ob_limbs[i], b->ob_limbs[i]);
        if (c != 0)
            return c;
    }
    return 0;
}

/* Python 3.x only understands rich comparisons; deccoeff_richcompare wraps
   _deccoeff_compare. */

static PyObject *
deccoeff_richcompare(PyObject *v, PyObject *w, int op)
{
    deccoeff *a, *b;
    PyObject *z = NULL;
    if (!compatible_with_deccoeff(v)) {
        Py_INCREF(Py_NotImplemented);
        z = Py_NotImplemented;
    }
    else if ((a = convert_to_deccoeff(v)) != NULL) {
        if (!compatible_with_deccoeff(w)) {
            Py_INCREF(Py_NotImplemented);
            z = Py_NotImplemented;
        }
        else if ((b = convert_to_deccoeff(w)) != NULL) {
            z = Py_CmpToRich(op, _deccoeff_compare(a, b));
            Py_DECREF(b);
        }
        Py_DECREF(a);
    }
    return z;
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
                        "Can't convert negative integer to " CLASS_NAME);
        return NULL;
    }
    z_size = scale_Py_ssize_t(a_size,
                              LOG2_10LQ * PyLong_SHIFT,
                              LOG2_10LP * LIMB_DIGITS);
    if (z_size == -1)
        PyErr_SetString(PyExc_OverflowError,
                        "Overflow in int to " CLASS_NAME " conversion\n");
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
    z_size = scale_Py_ssize_t(a_size,
                              LOG2_10UP * LIMB_DIGITS,
                              LOG2_10UQ * PyLong_SHIFT);
    if (z_size == -1)
        PyErr_SetString(PyExc_OverflowError,
                        "Overflow in " CLASS_NAME " to int conversion\n");
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

    /* not allowing subtypes */
    assert(type == &deccoeff_DeccoeffType);
    x = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O:" CLASS_NAME,
                                     kwlist, &x))
        return NULL;

    if (x == NULL)
        return (PyObject *)deccoeff_zero();
    else if (PyUnicode_Check(x)) {
        Py_UNICODE *s;
        Py_ssize_t s_len;
        s = PyUnicode_AS_UNICODE(x);
        s_len = PyUnicode_GET_SIZE(x);
        return (PyObject *)_deccoeff_from_unicode_and_size(s, s_len);
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
_deccoeff_length(deccoeff *v)
{
    Py_ssize_t v_size;
    v_size = Py_SIZE(v);
    if (v_size == 0)
        return 0;
    return limb_msd(v->ob_limbs[v_size-1]) + (v_size-1) * LIMB_DIGITS;
}

static PyObject *
deccoeff_digit_length(deccoeff *v)
{
    return PyLong_FromSsize_t(_deccoeff_length(v));
}

static void
deccoeff_dealloc(PyObject *v)
{
    Py_TYPE(v)->tp_free(v);
}

#if SIZEOF_LONG == 4
#define HASH_BITS 32
#define HASH_SHIFT 13
#define HASH_MASK ((1UL<<HASH_SHIFT) - 1)
#define HASH_START 1887730231UL
#elif SIZEOF_LONG == 8
#define HASH_BITS 64
#define HASH_SHIFT 25
#define HASH_MASK ((1UL<<HASH_SHIFT) - 1)
#define HASH_START 16569463434574008864UL
#else
#error "Expecting sizeof(long) to be 4 or 8"
#endif

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
        x ^= (unsigned long)limb_hash(*--v_top);
    }

    y = (long)x;
    return (y == -1) ? -2 : y;
}

DECCOEFF_WRAP_BINOP(deccoeff_add, _deccoeff_add)
DECCOEFF_WRAP_BINOP(deccoeff_subtract, _deccoeff_subtract)
DECCOEFF_WRAP_BINOP(deccoeff_multiply, _deccoeff_multiply)
DECCOEFF_WRAP_BINOP(deccoeff_remainder, _deccoeff_remainder)
DECCOEFF_WRAP_BINOP(deccoeff_divmod, _deccoeff_divmod)
DECCOEFF_WRAP_BINOP(deccoeff_floor_divide, _deccoeff_floor_divide)

static PyMethodDef deccoeff_methods[] = {
    {"digit_length", (PyCFunction)deccoeff_digit_length, METH_NOARGS,
     "Number of digits."},
    {NULL, NULL}
};

static PyMappingMethods deccoeff_as_mapping = {
    0,                                      /*mp_length*/
    (binaryfunc)deccoeff_subscript,         /*mp_subscript*/
    0, /*mp_ass_subscript*/
};

static PyNumberMethods deccoeff_as_number = {
    deccoeff_add,                           /*nb_add*/
    deccoeff_subtract,                      /*nb_subtract*/
    deccoeff_multiply,                      /*nb_multiply*/
    deccoeff_remainder,                     /*nb_remainder*/
    deccoeff_divmod,                        /*nb_divmod*/
    deccoeff_power,                         /*nb_power*/
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

    sz = _deccoeff_length(v);

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
            *--p = limb_to_digit(
                limb_rshift(&limb_value,
                            limb_value, 1, LIMB_ZERO));
    }
    /* most significant limb_t */
    limb_value = *limb_pointer;
    assert(limb_bool(limb_value));
    while (limb_bool(limb_value))
        *--p = limb_to_digit(
            limb_rshift(&limb_value, limb_value, 1, LIMB_ZERO));
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


/*****************
 * _Decimal type *
 *****************/

/* class name */
#define DECIMAL_NAME "_Decimal"

/* A finite decimal object needs a sign, a coefficient and an exponent.  An
   infinity has a sign and nothing more; the coefficient and exponent are
   ignored.  A (quiet or signalling) nan has a sign, and may carry additional
   information in the coefficient.  The exponent is not used. */

#define DEC_FLAGS_NEG (1<<0)
#define DEC_FLAGS_NONZERO (1<<1) /* currently unused, but may want this later */
#define DEC_FLAGS_SPECIAL (1<<2)
#define DEC_FLAGS_INF  (1<<3)
#define DEC_FLAGS_NAN (1<<4)
#define DEC_FLAGS_SNAN (1<<5)
#define DEC_FLAGS_QNAN (1<<6)

#define FINITE_FLAGS 0
#define INF_FLAGS (DEC_FLAGS_SPECIAL | DEC_FLAGS_INF)
#define QNAN_FLAGS (DEC_FLAGS_SPECIAL | DEC_FLAGS_NAN | DEC_FLAGS_QNAN)
#define SNAN_FLAGS (DEC_FLAGS_SPECIAL | DEC_FLAGS_NAN | DEC_FLAGS_SNAN)

/* note that FINITE_FLAGS < INF_FLAGS < SNAN_FLAGS < QNAN_FLAGS,
   corresponding to the standard total ordering of Decimals */

typedef int dec_flag_t;
typedef Py_ssize_t exp_t;

typedef struct {
    PyObject_HEAD
    deccoeff *dec_coeff;   /* coefficient, or NaN payload; NULL for infinity */
    PyLongObject *dec_exp; /* exponent: NULL for an infinity or NaN */
    dec_flag_t dec_flags;  /* flags describing sign and number class */
} _Decimal;

static PyTypeObject deccoeff__DecimalType;

/* create new _Decimal (or instance of a subclass of _Decimal); no
   validation or conversion---just allocate memory and fill the slots */

static _Decimal *
__Decimal_new(PyTypeObject *type, dec_flag_t flags, deccoeff *coeff,
              PyLongObject *exp)
{
    _Decimal *self;

    /* sanity checks */
    assert(
        ((flags | 1) == (FINITE_FLAGS | 1)) ||    /* finite (poss. zero) */
        ((flags | 1) == (INF_FLAGS | 1)) ||       /* infinity */
        ((flags | 1) == (QNAN_FLAGS | 1)) ||      /* quiet nan */
        ((flags | 1) == (SNAN_FLAGS | 1)));       /* signaling nan */

    /* incref coefficient and exponent, but only if they're non-NULL */
    if (flags & DEC_FLAGS_INF)
        assert(coeff == NULL);
    else {
        assert(coeff != NULL);
        Py_INCREF(coeff);
    }
    if (flags & DEC_FLAGS_SPECIAL)
        assert(exp == NULL);
    else {
        assert(exp != NULL);
        Py_INCREF(exp);
    }

    self = (_Decimal *)type->tp_alloc(type, 0);
    if (self == NULL)
        return NULL;
    self->dec_flags = flags;
    self->dec_coeff = coeff;
    self->dec_exp = exp;
    return self;
}

/* deallocate _Decimal instance */

static void
_Decimal_dealloc(PyObject *self)
{
    _Decimal *dec;
    dec = (_Decimal *)self;
    /* coefficient should be present iff decimal instance is not infinite */
    if (dec->dec_flags & DEC_FLAGS_INF)
        assert(dec->dec_coeff == NULL);
    else {
        assert(dec->dec_coeff != NULL);
        Py_DECREF(dec->dec_coeff);
    }
    /* exponent is present only for finite numbers */
    if (dec->dec_flags & DEC_FLAGS_SPECIAL)
        assert(dec->dec_exp == NULL);
    else {
        assert(dec->dec_exp != NULL);
        Py_DECREF(dec->dec_exp);
    }
    self->ob_type->tp_free(self);
}

/* Macros to create finite, infinite, qnan, and snan _Decimal instances */

#define FINITE_DECIMAL(type, sign, coeff, exp) (__Decimal_new(type, sign, coeff, exp))
#define INF_DECIMAL(type, sign) (__Decimal_new(type, sign | INF_FLAGS, NULL, NULL))
#define QNAN_DECIMAL(type, sign, payload) (__Decimal_new(type, sign | QNAN_FLAGS, payload, NULL))
#define SNAN_DECIMAL(type, sign, payload) (__Decimal_new(type, sign | SNAN_FLAGS, payload, NULL))

/* create a new finite _Decimal instance */

static PyObject *
_Decimal_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    /* create a new finite _Decimal instance, given a sequence consisting
       of its sign (an integer), coefficient and exponent */

    PyObject *ocoeff, *oexp;
    deccoeff *coeff;
    PyLongObject *exp;
    int sign;
    static char *kwlist[] = {"sign", "coeff", "exp", 0};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iOO:" DECIMAL_NAME, kwlist,
                                     &sign, &ocoeff, &oexp))
        return NULL;
    if (!(sign == 0 || sign == 1)) {
        PyErr_SetString(PyExc_ValueError, "sign should be 0 or 1");
        return NULL;
    }
    if (ocoeff->ob_type != &deccoeff_DeccoeffType) {
        PyErr_SetString(PyExc_TypeError, "coeff should have type " CLASS_NAME);
        return NULL;
    }
    coeff = (deccoeff *)ocoeff;
    if (oexp->ob_type != &PyLong_Type) {
        PyErr_SetString(PyExc_TypeError, "exp should have type int");
        return NULL;
    }
    exp = (PyLongObject *)oexp;
    return (PyObject *)FINITE_DECIMAL(type, sign, coeff, exp);
}

/* Create a _Decimal instance directly from a string; classmethod */

static PyObject *
_Decimal_from_str(PyTypeObject *cls, PyObject *arg)
{
    PyObject *result=NULL;
    deccoeff *coeff;
    Py_ssize_t ndigits;
    Py_UNICODE *coeff_start, *coeff_end, *s_end, *int_end, *exp_start, *s;
    int sign = 0, exp_sign = 0;
    PyObject *exp, *temp, *frac_digits;
    deccoeff *duexp;

    if (!PyUnicode_Check(arg)) {
        PyErr_SetString(PyExc_TypeError,
                        "Expected str instance");
        return NULL;
    }
    s = PyUnicode_AS_UNICODE(arg);
    s_end = s + PyUnicode_GET_SIZE(arg);

    /* Stage 1: parse and validate the string, identifying the points
       where the coefficients and exponent start and stop */

    /* optional sign */
    if (*s == '+')
        s++;
    else if (*s == '-') {
        s++;
        sign = 1;
    }

    switch(*s) {
    case 'n':
    case 'N':
        /* nan */
        s++;
        if (*s == 'a' || *s == 'A')
            s++;
        else
            goto parse_error;
        if (*s == 'n' || *s == 'N')
            s++;
        else
            goto parse_error;
        coeff_start = s;
        while ('0' <= *s && *s <= '9')
            s++;
        if (s != s_end)
            goto parse_error;

        coeff = _deccoeff_from_unicode_and_size(coeff_start, s-coeff_start);
        if (coeff != NULL) {
            result = (PyObject *)QNAN_DECIMAL(cls, sign, coeff);
            Py_DECREF(coeff);
        }
        break;
    case 's':
    case 'S':
        /* snan */
        s++;
        if (*s == 'n' || *s == 'N')
            s++;
        else
            goto parse_error;
        if (*s == 'a' || *s == 'A')
            s++;
        else
            goto parse_error;
        if (*s == 'n' || *s == 'N')
            s++;
        else
            goto parse_error;
        coeff_start = s;
        while ('0' <= *s && *s <= '9')
            s++;
        if (s != s_end)
            goto parse_error;

        coeff = _deccoeff_from_unicode_and_size(coeff_start, s-coeff_start);
        if (coeff != NULL) {
            result = (PyObject *)SNAN_DECIMAL(cls, sign, coeff);
            Py_DECREF(coeff);
        }
        break;
    case 'i':
    case 'I':
        /* inf[inity] */
        s++;
        if (*s == 'n' || *s == 'N')
            s++;
        else
            goto parse_error;
        if (*s == 'f' || *s == 'F')
            s++;
        else
            goto parse_error;
        if (*s == 'i' || *s == 'I') {
            s++;
            if (*s == 'n' || *s == 'N')
                s++;
            else
                goto parse_error;
            if (*s == 'i' || *s == 'I')
                s++;
            else
                goto parse_error;
            if (*s == 't' || *s == 'T')
                s++;
            else
                goto parse_error;
            if (*s == 'y' || *s == 'Y')
                s++;
            else
                goto parse_error;
        }
        /* end of string */
        if (s != s_end)
            goto parse_error;

        result = (PyObject *)INF_DECIMAL(cls, sign);
        break;
    default:
        /* numeric part: at least one digit, with an optional decimal point */
        coeff_start = s;
        while ('0' <= *s && *s <= '9')
            s++;
        int_end = s;
        if (*s == '.') {
            s++;
            while ('0' <= *s && *s <= '9')
                s++;
            coeff_end = s-1;
        }
        else
            coeff_end= s;

        ndigits = coeff_end - coeff_start;
        if (ndigits == 0)
            goto parse_error;

        /* [e <exponent>] */
        if (*s == 'e' || *s == 'E') {
            s++;
            if (*s == '+')
                s++;
            else if (*s == '-') {
                s++;
                exp_sign = 1;
            }
            exp_start = s;
            if (!('0' <= *s && *s <= '9'))
                goto parse_error;
            s++;
            while ('0' <= *s && *s <= '9')
                s++;
        }
        else
            exp_start = s;

        /* end of string */
        if (s != s_end)
            goto parse_error;

        /* parse exponent (without sign), returning a deccoeff */
        duexp = _deccoeff_from_unicode_and_size(exp_start, s-exp_start);
        if (duexp == NULL)
            return NULL;

        /* REF: duexp */

        /* Later we'll allow negative deccoeffs, and have the exponent
           be a deccoeff.   Later  :)
           For now we have to convert this to a Python long */
        exp = (PyObject *)deccoeff_long(duexp);
        Py_DECREF(duexp);
        if (exp == NULL)
            return NULL;

        /* REF: exp */

        /* adjust exponent:  include sign, and adjust by length
           of fractional part of input */
        if (exp_sign == 1) {
            /* negate exp */
            temp = PyNumber_Negative(exp);
            Py_DECREF(exp);
            exp = temp;
            if (exp == NULL)
                return NULL;
        }

        /* REF: exp */

        /* subtract frac_digits */
        frac_digits = PyLong_FromSize_t(coeff_end - int_end);
        if (frac_digits == NULL) {
            Py_DECREF(exp);
            return NULL;
        }

        /* REF: exp, frac_digits */

        temp = PyNumber_Subtract(exp, frac_digits);
        Py_DECREF(frac_digits);
        Py_DECREF(exp);
        exp = temp;
        if (exp == NULL)
            return NULL;

        /* REF: exp */

        /* get coefficient */
        coeff = _deccoeff_from_pointed_unicode_and_size(coeff_start,
                            coeff_end - coeff_start, int_end - coeff_start);

        if (coeff == NULL) {
            Py_DECREF(exp);
            return NULL;
        }

        result = (PyObject *)FINITE_DECIMAL(cls, sign, coeff,
                                              (PyLongObject *)exp);
        Py_DECREF(coeff);
        Py_DECREF(exp);
    }
    return result;

  parse_error:
    PyErr_SetString(PyExc_ValueError,
                    "invalid numeric string");
    return NULL;

}

/* Create a new finite _Decimal instance; classmethod */

static PyObject *
_Decimal_finite(PyTypeObject *cls, PyObject *args)
{
    PyObject *ocoeff, *oexp;
    deccoeff *coeff;
    PyLongObject *exp;
    int sign;

    if (!PyArg_ParseTuple(args, "iOO:" DECIMAL_NAME, &sign, &ocoeff, &oexp))
        return NULL;
    if (ocoeff->ob_type != &deccoeff_DeccoeffType) {
        PyErr_SetString(PyExc_TypeError, "coeff should have type " CLASS_NAME);
        return NULL;
    }
    coeff = (deccoeff *)ocoeff;
    if (oexp->ob_type != &PyLong_Type) {
        PyErr_SetString(PyExc_TypeError, "exp should have type int");
        return NULL;
    }
    exp = (PyLongObject *)oexp;
    if (!(sign == 0 || sign == 1)) {
        PyErr_SetString(PyExc_ValueError, "sign should be 0 or 1");
        return NULL;
    }
    return (PyObject *)FINITE_DECIMAL(cls, sign, coeff, exp);
}

/* Create a qNaN; classmethod */

static PyObject *
_Decimal_qNaN(PyTypeObject *cls, PyObject *args) {
    PyObject *opayload;
    deccoeff *payload;
    int sign;

    if (!PyArg_ParseTuple(args, "iO:" DECIMAL_NAME, &sign, &opayload))
        return NULL;
    if (opayload->ob_type != &deccoeff_DeccoeffType) {
        PyErr_SetString(PyExc_TypeError, "payload should have type " CLASS_NAME);
        return NULL;
    }
    payload = (deccoeff *)opayload;
    if (!(sign == 0 || sign == 1)) {
        PyErr_SetString(PyExc_ValueError, "sign should be 0 or 1");
        return NULL;
    }
    return (PyObject *)QNAN_DECIMAL(cls, sign, payload);
}

/* Create an sNaN; classmethod */

static PyObject *
_Decimal_sNaN(PyTypeObject *cls, PyObject *args) {
    PyObject *opayload;
    deccoeff *payload;
    int sign;

    if (!PyArg_ParseTuple(args, "iO:" DECIMAL_NAME, &sign, &opayload))
        return NULL;
    if (opayload->ob_type != &deccoeff_DeccoeffType) {
        PyErr_SetString(PyExc_TypeError, "payload should have type " CLASS_NAME);
        return NULL;
    }
    payload = (deccoeff *)opayload;
    if (!(sign == 0 || sign == 1)) {
        PyErr_SetString(PyExc_ValueError, "sign should be 0 or 1");
        return NULL;
    }
    return (PyObject *)SNAN_DECIMAL(cls, sign, payload);
}

/* Create an infinity; classmethod */

static PyObject *
_Decimal_inf(PyTypeObject *cls, PyObject *args) {
    int sign;

    if (!PyArg_ParseTuple(args, "i:" DECIMAL_NAME, &sign))
        return NULL;
    if (!(sign == 0 || sign == 1)) {
        PyErr_SetString(PyExc_ValueError, "sign should be 0 or 1");
        return NULL;
    }
    return (PyObject *)INF_DECIMAL(cls, sign);
}

/* Return the sign of any _Decimal instance */

static PyObject *
_Decimal_getsign(_Decimal *self, void *closure)
{
    long sign;
    sign = (long)(self->dec_flags) & DEC_FLAGS_NEG;
    return PyLong_FromLong(sign);
}

/* Return the coefficient of a finite _Decimal */

static PyObject *
_Decimal_getcoeff(_Decimal *self, void *closure)
{
    if (self->dec_flags & DEC_FLAGS_SPECIAL) {
        PyErr_SetString(PyExc_ValueError,
                        "infinity or NaN has no coefficient");
        return NULL;
    }
    Py_INCREF(self->dec_coeff);
    return (PyObject *)self->dec_coeff;
}

/* Return the payload of a NaN */

static PyObject *
_Decimal_getpayload(_Decimal *self, void *closure)
{
    if (self->dec_flags & DEC_FLAGS_NAN) {
        Py_INCREF(self->dec_coeff);
        return (PyObject *)self->dec_coeff;
    }
    else {
        PyErr_SetString(PyExc_ValueError,
                        "argument is not a NaN");
        return NULL;
    }
}

/* Return the exponent of a finite _Decimal instance */

static PyObject *
_Decimal_getexp(_Decimal *self, void *closure)
{
    if (self->dec_flags & DEC_FLAGS_SPECIAL) {
        PyErr_SetString(PyExc_ValueError,
                        "infinity or NaN has no exponent");
        return NULL;
    }
    Py_INCREF(self->dec_exp);
    return (PyObject *)(self->dec_exp);
}

/* Return True if the given Decimal is special (an infinity or NaN),
   False otherwise. */

static PyObject *
_Decimal_getspecial(_Decimal *self, void *closure)
{
    if (self->dec_flags & DEC_FLAGS_SPECIAL)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}


static PyObject *
_Decimal_is_finite(_Decimal *self)
{
    if (self->dec_flags & DEC_FLAGS_SPECIAL)
        Py_RETURN_FALSE;
    else
        Py_RETURN_TRUE;
}

static PyObject *
_Decimal_is_infinite(_Decimal *self)
{
    if (self->dec_flags & DEC_FLAGS_INF)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
_Decimal_is_nan(_Decimal *self)
{
    if (self->dec_flags & DEC_FLAGS_NAN)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
_Decimal_is_qnan(_Decimal *self)
{
    if (self->dec_flags & DEC_FLAGS_QNAN)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
_Decimal_is_snan(_Decimal *self)
{
    if (self->dec_flags & DEC_FLAGS_SNAN)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
_Decimal_is_signed(_Decimal *self)
{
    if (self->dec_flags & DEC_FLAGS_NEG)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
_Decimal_is_zero(_Decimal *self)
{
    if ((self->dec_flags & DEC_FLAGS_SPECIAL) ||
        deccoeff_bool(self->dec_coeff))
        Py_RETURN_FALSE;
    else
        Py_RETURN_TRUE;
}

static _Decimal *
_Decimal_copy(_Decimal *self)
{
    return __Decimal_new(&deccoeff__DecimalType, self->dec_flags,
                         self->dec_coeff, self->dec_exp);
}

static _Decimal *
_Decimal_copy_negate(_Decimal *self)
{
    return __Decimal_new(&deccoeff__DecimalType, self->dec_flags ^ 1,
                         self->dec_coeff, self->dec_exp);
}

static _Decimal *
_Decimal_copy_abs(_Decimal *self)
{
    return __Decimal_new(&deccoeff__DecimalType, self->dec_flags & (~1),
                         self->dec_coeff, self->dec_exp);
}

static PyMethodDef _Decimal_methods[] = {
    {"_finite", (PyCFunction)_Decimal_finite, METH_VARARGS|METH_CLASS, " "},
    {"_qnan", (PyCFunction)_Decimal_qNaN, METH_VARARGS|METH_CLASS, " "},
    {"_snan", (PyCFunction)_Decimal_sNaN, METH_VARARGS|METH_CLASS, " "},
    {"_inf", (PyCFunction)_Decimal_inf, METH_VARARGS|METH_CLASS, " "},
    {"from_str", (PyCFunction)_Decimal_from_str, METH_O|METH_CLASS, " "},
    {"is_finite", (PyCFunction)_Decimal_is_finite, METH_NOARGS, " "},
    {"is_infinite", (PyCFunction)_Decimal_is_infinite, METH_NOARGS, " "},
    {"is_nan", (PyCFunction)_Decimal_is_nan, METH_NOARGS, " "},
    {"is_qnan", (PyCFunction)_Decimal_is_qnan, METH_NOARGS, " "},
    {"is_snan", (PyCFunction)_Decimal_is_snan, METH_NOARGS, " "},
    {"is_signed", (PyCFunction)_Decimal_is_signed, METH_NOARGS, " "},
    {"is_zero", (PyCFunction)_Decimal_is_zero, METH_NOARGS, " "},
    {"copy", (PyCFunction)_Decimal_copy, METH_NOARGS, " "},
    {"copy_negate", (PyCFunction)_Decimal_copy_negate, METH_NOARGS, " "},
    {"copy_abs", (PyCFunction)_Decimal_copy_abs, METH_NOARGS, " "},
    {NULL, NULL}
};

static PyGetSetDef _Decimal_getsetters[] = {
    {"_sign", (getter)_Decimal_getsign, NULL, "sign", NULL},
    {"_int", (getter)_Decimal_getcoeff, NULL,
     "coefficient (invalid for NaNs and infinites)", NULL},
    {"_exp", (getter)_Decimal_getexp, NULL,
     "exponent (invalid for NaNs and infinities)", NULL},
    {"_payload", (getter)_Decimal_getpayload, NULL,
     "payload of a NaN (invalid for non-NaNs)", NULL},
    {"_is_special", (getter)_Decimal_getspecial, NULL,
     "True for infinities and NaNs, false otherwise", NULL},
    {NULL}
};

static PyTypeObject deccoeff__DecimalType = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    MODULE_NAME "." DECIMAL_NAME,             /* tp_name */
    sizeof(_Decimal),                       /* tp_basicsize */
    0,                          /* tp_itemsize */
    _Decimal_dealloc,                       /* tp_dealloc */
    0, /* tp_print */
    0, /* tp_getattr */
    0, /* tp_setattr */
    0, /* tp_compare */
    0,                /* tp_repr */
    0,                    /* tp_as_number */
    0, /* tp_as_sequence */
    0,                   /* tp_as_mapping */
    0,                /* tp_hash */
    0, /* tp_call */
    0,                 /* tp_str */
    0, /* tp_getattro */
    0, /* tp_setattro */
    0, /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,                     /* tp_flags */
    "support for Decimal type",                     /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0,      /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    _Decimal_methods,                       /* tp_methods */
    0, /* tp_members */
    _Decimal_getsetters, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    0, /* tp_init */
    0, /* tp_alloc */
    _Decimal_new,                           /* tp_new */
    PyObject_Del,                           /* tp_free */
};

static PyMethodDef deccoeff_module_methods[] = {
    {NULL, NULL}
};

static struct PyModuleDef deccoeff_module = {
    PyModuleDef_HEAD_INIT,
    "deccoeff",
    "class for decimal integer arithmetic; support for decimal module",
    -1,
    deccoeff_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_deccoeff(void)
{
    PyObject *m;
    int check, i;
    limb_t power_of_ten;

    /* initialize powers_of_ten array */
    power_of_ten = 1;
    for (i = 0; i < LIMB_DIGITS; i++) {
        powers_of_ten[i] = power_of_ten;
        power_of_ten *= 10;
    }

    if (PyType_Ready(&deccoeff_DeccoeffType) < 0)
        return NULL;

    if (PyType_Ready(&deccoeff__DecimalType) < 0)
        return NULL;

    m = PyModule_Create(&deccoeff_module);
    if (m == NULL)
        return NULL;

    Py_INCREF(&deccoeff_DeccoeffType);
    check = PyModule_AddObject(m, CLASS_NAME,
                               (PyObject *) &deccoeff_DeccoeffType);
    if (check == -1)
        return NULL;

    Py_INCREF(&deccoeff__DecimalType);
    check = PyModule_AddObject(m, DECIMAL_NAME,
                               (PyObject *) &deccoeff__DecimalType);
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
