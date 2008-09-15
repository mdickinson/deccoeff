/*
 *  To do
 *  -----
 *  Replace limbs.c/limbs.h with an implementation involving strings (or
 *    structs), just to check that _decimal.c is properly independent of the
 *    exact implementation.
 *  Karatsuba multiplication
 *  docstrings
 *  make sure that mixed-type arithmetic raises a suitable exception
 *  fix ob_limbs[0]; it's not valid C89
 *  three-argument pow
 *  make slicing quicker if the shift is a multiple of LIMB_DIGITS
 *  check for overflow in slicing
 *  consider semantics for boolean operations: one possibility would
 *      be to consider any nonzero digit to be equivalent to 1
 *  implement powering algorithm based on base 10
 *  write Deccoeff-specific tests
 *  export LIMB_DIGITS to Python
 *  make sure single-limb adds, subs, muls are as fast as possible
 *
 */

/*
 * Module implementing high-precision natural number arithmetic in a
 * decimal base.  As the name 'deccoeff' suggests, these numbers are
 * intended to be used as the coefficients for Decimal instances.
 *
 * Author: Mark Dickinson.  Licensed to PSF under a Contributor
 * Agreement.
 */


/* Various notes:

   Right-hand operands to shifts are assumed to be binary.
   Similarly, indices and slices are assumed to be given in binary.
   This make sense assuming that the decimal exponent is going to
   be stored in binary.

*/
#include "Python.h"
#include "longintrepr.h"
#include "limb.h"
#include "string.h"

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

typedef limb_t *limbs;
typedef const limb_t *const_limbs;

static bool
limbs_inc(limbs res, const_limbs a, Py_ssize_t n)
{
	Py_ssize_t i;
	bool carry;
	carry = true;
	for (i=0; i < n; i++)
		carry = limb_incc(res+i, a[i], carry);
	return carry;
}

static bool
limbs_incc(limbs res, const_limbs a, Py_ssize_t n, bool carry)
{
	Py_ssize_t i;
	for (i=0; i < n; i++)
		carry = limb_incc(res+i, a[i], carry);
	return carry;
}

static bool
limbs_dec(limbs res, const_limbs a, Py_ssize_t n)
{
	Py_ssize_t i;
	bool carry;
	carry = true;
	for (i=0; i < n; i++)
		carry = limb_sbb(res+i, a[i], LIMB_ZERO, carry);
	return carry;
}

static bool
limbs_decc(limbs res, const_limbs a, Py_ssize_t n, bool carry)
{
	Py_ssize_t i;
	for (i=0; i < n; i++)
		carry = limb_sbb(res+i, a[i], LIMB_ZERO, carry);
	return carry;
}

static void
limbs_copy(limbs res, const_limbs a, Py_ssize_t n)
{
	Py_ssize_t i;
	for (i=0; i < n; i++)
		res[i] = a[i];
}

/* add n-limb numbers a and b, producing an n-limb result res and a carry */

static bool
limbs_add(limbs res, const_limbs a, const_limbs b, Py_ssize_t n)
{
	Py_ssize_t i;
	bool carry;
	carry = false;
	for (i=0; i < n; i++)
		carry = limb_adc(res+i, a[i], b[i], carry);
	return carry;
}

/* a-b, a and b n-limb numbers.  n-limb result in res; return carry. */

static bool
limbs_sub(limbs res, const_limbs a, const_limbs b, Py_ssize_t n)
{
	Py_ssize_t i;
	bool carry;
	carry = false;
	for (i=0; i < n; i++)
		carry = limb_sbb(res+i, a[i], b[i], carry);
	return carry;
}

/* multiply m-limb number a by single limb x, getting m-digit result res and
   an extra high limb. */

static limb_t
limbs_mul1(limbs res, const_limbs a, Py_ssize_t n, limb_t x)
{
	limb_t high;
	Py_ssize_t i;
	high = LIMB_ZERO;
	for (i=0; i < n; i++)
		high = limb_fmaa(res+i, a[i], x, high, LIMB_ZERO);
	return high;
}

/* multiply m-digit a by n-digit b, getting m+n-digit result res.
   res should not overlap either of the inputs. */

static void
limbs_mul(limbs res, const_limbs a, Py_ssize_t m, const_limbs b, Py_ssize_t n)
{
	Py_ssize_t i, j;
	limb_t hiword;
	for (j=0; j < n; j++)
		res[j] = LIMB_ZERO;
	for (i=0; i < m; i++) {
		hiword = LIMB_ZERO;
		for (j=0; j < n; j++)
			hiword = limb_fmaa(res+i+j, a[i], b[j],
					   res[i+j], hiword);
		res[i+j] = hiword;
	}
}

/* divide m-limb number a by single limb x, giving m-digit quotient res and
   returning the remainder */

static limb_t
limbs_div1(limbs res, const_limbs a, Py_ssize_t m, limb_t x)
{
	limb_t high;
	Py_ssize_t i;
	high = LIMB_ZERO;
	for (i = m-1; i >= 0; i--)
		res[i] = limb_div(&high, high, a[i], x);
	return high;
}

/* print limbs to stdout, most significant first; used for debugging only */

static void
limbs_print(const_limbs a, Py_ssize_t n)
{
	Py_ssize_t i;
	printf("[");
	for (i = n-1; i >= 0; i--) {
		printf("%0*d", LIMB_DIGITS, a[i]);
		if (i != 0)
			printf(" ");
	}
	printf("]");
}

/* divide m-limb a by n-limb b, giving an (m-n+1)-limb quotient and n-limb
   remainder.  Assumes that the top digit of b is nonzero.  w provides m+n+1
   digits of workspace. */

static void
limbs_div(limbs quot, limbs rem, const_limbs a, Py_ssize_t m, const_limbs b,
	  Py_ssize_t n, limbs w)
{
	limb_t scale, top, a_top, b_top, q, dummy;
	limbs aa, bb;
	bool carry;
	Py_ssize_t j;

	/* top limb of b should be nonzero; a should contain at least as many
	   limbs as b */
	assert(m >= n && n > 0 && !limb_eq(b[n-1], LIMB_ZERO));

	/* compute scale factor for normalization: floor(LIMB_BASE /
	   (b_top+1)) */
	carry = limb_inc(&scale, b[n-1]);
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
	if (top == 0 && aa[m-1] < bb[n-1])
		quot[m-n] = 0;
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
		carry = limb_incc(&top, top, carry);
		assert(!carry);
		assert(limb_le(a_top, top));
		/* correct if necessary */
		while (limb_lt(a_top, top)) {
			carry = limbs_add(aa, aa, bb, n);
			carry = limb_incc(&a_top, a_top, carry);
			assert(!carry);
			carry = limb_dec(&q, q);
			assert(!carry);
		}
		quot[j] = q;
	}
	/* at the end of this loop, the bottom n limbs of a give the
	   remainder, which needs to be divided through by the scale factor.
	   XXX this step can be skipped when all we want is the quotient. */
	top = limbs_div1(rem, aa, n, scale);
	assert(limb_eq(top, LIMB_ZERO));
}

/* shift m-limb number right n digits (not limbs!), 0 <= n <= LIMB_DIGITS*m */

static void
limbs_rshift(limbs res, const_limbs a, Py_ssize_t m, Py_ssize_t n)
{
	Py_ssize_t n_limbs, n_digits, i;
	limb_t limb_top, limb_bot, rem;

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
		res[i] = rem + limb_bot;
		rem = limb_top;
	}
	res[i] = rem;
}

/* res = a[m:n], 0 <= m < n <= LIMB_DIGITS*a_size. res should have length
   ceiling((n-m)/LIMB_DIGITS). */

/* number of limbs of a needed is (n-1)/LIMB_DIGITS - m/LIMB_DIGITS + 1
   number of limbs in result is (n-1-m)/LIMB_DIGITS + 1.  These two
   numbers are the same iff (n-1) % LIMB_DIGITS >= m % LIMB_DIGITS; otherwise
   we need to use an extra limb of a.

   Ex: a[0:9] needs only one limb of a
       a[1:10] needs two limbs of a, one limb result
       a[2:11] needs two limbs of a, one limb result

       a[0:8] needs only one limb of a
       a[1:9] needs only one limb of a
       a[2:10] needs two limbs of a

       a[0:10] needs two limbs
       a[1:11] needs two limbs

       a[0:2] one limb
       a[7:9] one limb
       a[8:10] two
       a[9:11] one

 */

static void
limbs_slice(limbs res, const_limbs a, Py_ssize_t m, Py_ssize_t n)
{
	Py_ssize_t m_limbs, m_digits, res_limbs, res_digits, i;
	limb_t out, limb_bot, limb_top;

	/* number of limbs of result is (n-1-m)/LIMB_DIGITS + 1 */
	/* want to know whether (n-1-m) % LIMB_DIGITS > LIMB_DIGITS - m % LIMB_DIGITS */

	/* if a % LIMB_DIGITS < b % LIMB_DIGITS then a-b produces a carry; and
	   (a-b) % LIMB_DIGITS is at least LIMB_DIGITS - b % LIMB_DIGITS.  if
	   a % LIMB_DIGITS >= b % LIMB_DIGITS then a-b doesn't carry, and
	   (a-b) % LIMB_DIGITS is strictly less than LIMB_DIGITS - b %
	   LIMB_DIGITS.
	*/

	m_limbs = m / LIMB_DIGITS;
	m_digits = m % LIMB_DIGITS;
	res_limbs = (n-1-m) / LIMB_DIGITS;
	res_digits = (n-1-m) % LIMB_DIGITS + 1;
	out = limb_sar(a[m_limbs++], m_digits);
	for (i = 0; i < res_limbs; i++) {
		limb_top = limb_split(&limb_bot, a[m_limbs++], m_digits);
		res[i] = out + limb_bot;
		out = limb_top;
	}
	if (res_digits > LIMB_DIGITS - m_digits) {
		limb_bot = limb_sal(a[m_limbs++], LIMB_DIGITS - m_digits);
		out += limb_bot;
	}
	res[i] = limb_low(out, res_digits);
}

/* shift m-limb number left n digits (shifting zeros in); assumes that
   res has size m + ceiling(n/LIMB_DIGITS) */

static void
limbs_lshift(limbs res, const_limbs a, Py_ssize_t m, Py_ssize_t n)
{
	Py_ssize_t n_limbs, n_digits, i;
	limb_t limb_top, tmp;

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

	/* shift bottom digit up;  keep low part */
	assert(i == n_limbs);
	res[i++] = limb_splitl(&limb_top, a[0], n_digits);
	for (; i < n_limbs+m; i++) {
		res[i] = limb_top + limb_splitl(&tmp, a[i-n_limbs],
						n_digits);
		limb_top = tmp;
	}
	res[i] = limb_top;
}

/**************************
 * deccoeff : definitions *
 **************************/

typedef struct {
  PyObject_VAR_HEAD
  limb_t ob_limbs[0];
} deccoeff;

/* We place an upper bound MAX_DIGITS on the number of decimal digits (*not*
   the number of limbs).  MAX_DIGITS should fit into a Py_ssize_t, so that
   length and indexing always make sense. */

/* XXX add code to make sure that MAX_DIGITS fits in a Py_ssize_t.
   Or possibly to make sure that 2*MAX_DIGITS fits in a Py_ssize_t. */

#define MAX_DIGITS 1000000000

static deccoeff *deccoeff_const_zero;
static deccoeff *deccoeff_const_one;
static deccoeff *deccoeff_PyLong_BASE;

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

/* return nonzero value iff v has at most MAX_DIGITS digits. */

static int
deccoeff_checksize(deccoeff *v)
{
	Py_ssize_t v_size;
	v_size = Py_SIZE(v);
	return (v_size < (MAX_DIGITS-1)/LIMB_DIGITS+1 ||
		(v_size == (MAX_DIGITS-1)/LIMB_DIGITS+1 &&
		 limb_dsr(v->ob_limbs[v_size-1]) <=
		 (MAX_DIGITS-1)%LIMB_DIGITS+1));
}

/* return a new reference to the zero deccoeff */

static deccoeff *
deccoeff_zero(void)
{
	Py_INCREF(deccoeff_const_zero);
	return deccoeff_const_zero;
}

/* create a deccoeff from a C unsigned long. */

static deccoeff *
deccoeff_from_ulong(unsigned long x)
{
	unsigned long xcopy;
	Py_ssize_t z_size, i;
	deccoeff *z;

	/* wasteful method of figuring out number of limbs needed */
	z_size = 0;
	xcopy = x;
	while(xcopy > 0) {
		limb_from_ulong(&xcopy, xcopy);
		z_size++;
	}
	/* allocate space */
	z = _deccoeff_new(z_size);
	if (z==NULL)
		return NULL;
	/* compute result */
	for (i=0; i < z_size; i++)
		z->ob_limbs[i] = limb_from_ulong(&x, x);
	return z;
}

/***************************
 * Arithmetic on deccoeffs *
 ***************************/

/* addition: raise OverflowError if result has more than MAX_DIGITS digits */

static deccoeff *
deccoeff_add(deccoeff *a, deccoeff *b)
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
	z = deccoeff_normalize(z);
	if (deccoeff_checksize(z))
		return z;
	Py_DECREF(z);
	PyErr_SetString(PyExc_OverflowError,
			"Too many digits in sum");
	return NULL;
}

/* subtraction: return a - b; raise OverflowError if the result is negative */

static deccoeff *
deccoeff_subtract(deccoeff *a, deccoeff *b)
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

/* multiplication: return a * b.  Raise OverflowError if result too large. */

static deccoeff *
deccoeff_multiply(deccoeff *a, deccoeff *b)
{
	Py_ssize_t a_size, b_size;
	deccoeff *z;
	a_size = Py_SIZE(a);
	b_size = Py_SIZE(b);
	z = _deccoeff_new(a_size + b_size);
	if (z == NULL)
		return NULL;
	limbs_mul(z->ob_limbs, a->ob_limbs, a_size, b->ob_limbs, b_size);
	z = deccoeff_normalize(z);
	if (deccoeff_checksize(z))
		return z;
	Py_DECREF(z);
	PyErr_SetString(PyExc_OverflowError, "product has too many digits");
	return NULL;
}

/* division of a by b: returns quotient and puts remainder in *r.  On
   failure, both the quotient and the remainder are NULL.  Raises
   ZeroDivisionError on division by zero. */

static deccoeff *
_deccoeff_divmod(deccoeff **r, deccoeff *a, deccoeff *b) {
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
deccoeff_remainder(deccoeff *a, deccoeff *b) {
	deccoeff *quot, *rem;
	quot = _deccoeff_divmod(&rem, a, b);
	if (rem == NULL)
		return NULL;
	Py_DECREF(quot);
	return rem;
}

/* divmod: raises ZeroDivisionError if b is zero */

static PyObject *
deccoeff_divmod(deccoeff *a, deccoeff *b) {
	deccoeff *quot, *rem;

	quot = _deccoeff_divmod(&rem, a, b);
	if (quot == NULL)
		return NULL;
	return Py_BuildValue("OO", (PyObject *)quot, (PyObject *)rem);
}

/* negation: succeeds only for a == 0 */

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

/* left shift; second operand should be a nonnegative Python integer, not a
   deccoeff. raises ValueError if second operand is negative, and
   OverflowError if result has more than MAX_DIGITS digits. */

static deccoeff *
deccoeff_lshift(deccoeff *a, PyObject *b) {
	Py_ssize_t n, a_size;
	deccoeff *z;
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
		return deccoeff_zero();
	if (n >= MAX_DIGITS)
		goto Overflow;
	z = _deccoeff_new(a_size + (n+LIMB_DIGITS-1) / LIMB_DIGITS);
	if (z==NULL)
		return NULL;
	limbs_lshift(z->ob_limbs, a->ob_limbs, a_size, n);
	z = deccoeff_normalize(z);
	if (deccoeff_checksize(z))
		return z;
	Py_DECREF(z);
  Overflow:
	PyErr_SetString(PyExc_OverflowError,
			"too many digits in result");
	return NULL;
}

/* right shift; second operand is a Python integer, not a deccoeff.  Raises
   ValueError if second operand is negative. */

static deccoeff *
deccoeff_rshift(deccoeff *a, PyObject *b) {
	Py_ssize_t n, a_size, shift;
	deccoeff *z;
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
		return deccoeff_zero();
	z = _deccoeff_new(a_size - shift);
	if (z==NULL)
		return NULL;
	limbs_rshift(z->ob_limbs, a->ob_limbs, a_size, n);
	return deccoeff_normalize(z);
}

/* slice */

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
					"step not supported in deccoeff slice");
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

/* floor division */

static deccoeff *
deccoeff_floor_divide(deccoeff *a, deccoeff *b) {
	deccoeff *quot, *rem;
	quot = _deccoeff_divmod(&rem, a, b);
	if (quot == NULL)
		return NULL;
	Py_DECREF(rem);
	return quot;
}

/* compare: return -1 if a < b, 0 if a == b, 1 if a > b.  Always succeeds. */

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

/* Create a deccoeff from a Python long integer */

/* Naive algorithm; needs optimization. But not before testing with a range of
   limb_t sizes. */

static deccoeff *
deccoeff_from_PyLong(deccoeff *self, PyObject *o)
{
	Py_ssize_t a_size;
	PyLongObject *a;
	digit *a_start, *a_end;
	deccoeff *result, *result2, *addend;

	if (!PyLong_Check(o)) {
		PyErr_SetString(PyExc_TypeError,
				"argument must be a long");
		return NULL;
	}

	a = (PyLongObject *)o;
	a_size = Py_SIZE(a);

	if (a_size < 0) {
		PyErr_SetString(PyExc_OverflowError,
				"Can't convert negative integer to a deccoeff");
		return NULL;
	}

	result = deccoeff_zero();

	a_start = a->ob_digit;
	a_end = a_start+a_size;
	while (a_end > a_start) {

		/* multiply result by PyLong_BASE, freeing old reference */
		result2 = (deccoeff *)deccoeff_multiply(deccoeff_PyLong_BASE, result);
		Py_DECREF(result);
		result = result2;
		if (result == NULL)
			return NULL;

		/* add PyLong digit to result */
		addend = deccoeff_from_ulong((unsigned long)(*--a_end));
		if (addend == NULL) {
			Py_DECREF(result);
			return NULL;
		}

		result2 = (deccoeff *)deccoeff_add(result, addend);
		Py_DECREF(result);
		Py_DECREF(addend);
		result = result2;
		if (result == NULL)
			return NULL;
	}
	return result;
}

/* we also need (eventually) a routine to convert a deccoeff back into a
   Python long.  I'm resisting writing this for the moment, since the lack of
   it provides incentive to keep all computations in decimal. */

static PyObject *
deccoeff_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	const char* s;
	char *start, *end;
	int lens;
	int nlimbs, digits_in_limb;
	limb_t accum, *limb_pointer;
	char c;
	deccoeff *new;

	if (!PyArg_ParseTuple(args, "s#:" CLASS_NAME, &s, &lens))
		return NULL;

	if (lens > MAX_DIGITS) {
		PyErr_SetString(PyExc_OverflowError,
				"too many digits");
		return NULL;
	}

	nlimbs = (lens+LIMB_DIGITS-1)/LIMB_DIGITS;
	digits_in_limb = (lens+LIMB_DIGITS-1) % LIMB_DIGITS + 1;

	new = _deccoeff_new(nlimbs);
	if (new == NULL)
		return NULL;

	/* iterate through digits, from most significant to least */
	start = (char *)s;
	end = (char *)s + lens;
	accum = LIMB_ZERO;
	limb_pointer = new -> ob_limbs + nlimbs - 1;
	while (start < end) {
		c = *start++;
		if (!(c >= '0') || !(c <= '9')) goto Error;
		accum = limb_sal(accum, 1); /* shift left 1 place */
		accum = limb_setdigit(accum, 0, c);
		digits_in_limb--;
		if (digits_in_limb == 0) {
			digits_in_limb = LIMB_DIGITS;
			*limb_pointer-- = accum;
			accum = LIMB_ZERO;
		}
	}
	return (PyObject *)(deccoeff_normalize(new));

  Error:
	PyErr_Format(PyExc_ValueError,
		     "invalid literal for deccoeff(): %s", s);
	Py_DECREF(new);
	return NULL;
}

/* number of digits of a deccoeff.  Raises ValueError if the deccoeff is zero. */

static Py_ssize_t
deccoeff_length(deccoeff *v)
{
	Py_ssize_t v_size;

	v_size = Py_SIZE(v);
	/* for safety, we raise a ValueError for 0; there are too many
	   situations where returning 0 might turn out to be wrong; 0 almost
	   always needs special treatment anyway.

	   Another way to think of this function is as 1 more than the floor
	   of log_10(v); this again is undefined for v = 0, justifying the
	   ValueError. */
	if (v_size == 0) {
		PyErr_SetString(PyExc_ValueError,
				"refusing to find length of 0");
		return -1;
	}
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
		for (i=0; i < LIMB_DIGITS; i++) {
			*--p = limb_getdigit(limb_value, 0);
			limb_value = limb_sar(limb_value, 1);
		}
	}
	/* most significant limb_t */
	limb_value = *limb_pointer;
	assert(!limb_eq(limb_value, LIMB_ZERO));
	while (!limb_eq(limb_value, LIMB_ZERO)) {
		*--p = limb_getdigit(limb_value, 0);
		limb_value = limb_sar(limb_value, 1);
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
	{"from_int", (PyCFunction)deccoeff_from_PyLong, METH_O | METH_STATIC, "create from an integer"},
	{NULL, NULL}
};

static PyMappingMethods deccoeff_as_mapping = {
	(lenfunc)deccoeff_length,               /*mp_length*/
	(binaryfunc)deccoeff_subscript,         /*mp_subscript*/
	0, /*mp_ass_subscript*/
};

static PyNumberMethods deccoeff_as_number = {
	(binaryfunc) deccoeff_add,              /*nb_add*/
	(binaryfunc) deccoeff_subtract,         /*nb_subtract*/
	(binaryfunc) deccoeff_multiply,         /*nb_multiply*/
	(binaryfunc) deccoeff_remainder,        /*nb_remainder*/
	(binaryfunc) deccoeff_divmod,           /*nb_divmod*/
	0, /*nb_power*/
	(unaryfunc) deccoeff_negative,          /*nb_negative*/
	(unaryfunc) deccoeff_positive,          /*nb_positive*/
	(unaryfunc) deccoeff_positive,          /*nb_absolute*/
	(inquiry) deccoeff_bool,                /*nb_bool*/
	0, /*nb_invert*/
	(binaryfunc) deccoeff_lshift,           /*nb_lshift*/
	(binaryfunc) deccoeff_rshift,           /*nb_rshift*/
	0, /*nb_and*/
	0, /*nb_xor*/
	0, /*nb_or*/
	0, /*nb_int*/
	0, /*nb_long*/
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
	(binaryfunc)deccoeff_floor_divide,      /*nb_floor_divide*/
	0, /*nb_true_divide*/
	0, /*nb_inplace_floor_divide*/
	0, /*nb_inplace_true_divide*/
	0, /*nb_index*/
};

static PyTypeObject deccoeff_DeccoeffType = {
	PyVarObject_HEAD_INIT(&PyType_Type, 0)
	MODULE_NAME "." CLASS_NAME, /* tp_name */
	sizeof(deccoeff), /* tp_basicsize */
	sizeof(limb_t),            /* tp_itemsize */
	deccoeff_dealloc,                            /* tp_dealloc */
	0,                           /* tp_print */
	0,                           /* tp_getattr */
	0,                           /* tp_setattr */
	0,                           /* tp_compare */
	(reprfunc)deccoeff_repr,      /* tp_repr */
	&deccoeff_as_number,          /* tp_as_number */
	0,                           /* tp_as_sequence */
	&deccoeff_as_mapping,                           /* tp_as_mapping */
	(hashfunc)deccoeff_hash,      /* tp_hash */
	0,                           /* tp_call */
	(reprfunc)deccoeff_str,       /* tp_str */
	0,     /* tp_getattro */
	0,                           /* tp_setattro */
	0,                           /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT,   /* tp_flags */
	"Decimal integers",                 /* tp_doc */
	0,                           /* tp_traverse */
	0,                           /* tp_clear */
	(richcmpfunc)deccoeff_richcompare, /* tp_richcompare */
	0,                           /* tp_weaklistoffset */
	0,                           /* tp_iter */
	0,                           /* tp_iternext */
	deccoeff_methods,             /* tp_methods */
	0,            /* tp_members */
	0,             /* tp_getset */
	0,                           /* tp_base */
	0,                           /* tp_dict */
	0,                           /* tp_descr_get */
	0,                           /* tp_descr_set */
	0,                           /* tp_dictoffset */
	0,                           /* tp_init */
	0,         /* tp_alloc */
	deccoeff_new,                 /* tp_new */
	PyObject_Del,                /* tp_free */
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

#define ADD_CONST(m, name)                                  \
    if (PyModule_AddIntConstant(m, #name, name) < 0) return m

PyMODINIT_FUNC
PyInit__decimal(void)
{
	PyObject *m = NULL;    /* a module object */

	if (PyType_Ready(&deccoeff_DeccoeffType) < 0)
		return m;

	m = PyModule_Create(&_decimalmodule);
	if (m == NULL)
		return m;

	Py_INCREF(&deccoeff_DeccoeffType);
	PyModule_AddObject(m, CLASS_NAME, (PyObject *) &deccoeff_DeccoeffType);

	deccoeff_const_zero = deccoeff_from_ulong((unsigned long)0);
	if (deccoeff_const_zero == NULL)
		return NULL;

	deccoeff_const_one = deccoeff_from_ulong((unsigned long)1);
	if (deccoeff_const_one == NULL)
		return NULL;

	deccoeff_PyLong_BASE = deccoeff_from_ulong((unsigned long)PyLong_BASE);
	if (deccoeff_PyLong_BASE == NULL)
		return NULL;

	/* XXX do we need to deallocate deccoeff_zero and deccoeff_PyLong_BASE
	   somewhere?   Do multiple imports leak references? */
	return m;
}
