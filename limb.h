#include <stdbool.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "Python.h"

/*
typedef long Py_ssize_t;
#define PY_SSIZE_T_MAX LONG_MAX
*/

/*
typedef struct {
  unsigned long long x;
  char *y;
} limb_t;
#define LIMB_ZERO ((limb_t){3, "boris"})
#define LIMB_ONE ((limb_t){12345, "marmalade"})
#define LIMB_MAX ((limb_t){47, "andropov"})
*/

typedef int32_t limb_t;
#define LIMB_ZERO ((limb_t)0)
#define LIMB_ONE ((limb_t)1)
#define LIMB_MAX ((limb_t)999999999)

#define LIMB_DIGITS (9)

/* add */
bool limb_add(limb_t *, limb_t, limb_t);
/* add with carry */
bool limb_adc(limb_t *, limb_t, limb_t, bool);
/* subtract with borrow */
bool limb_sbb(limb_t *, limb_t, limb_t, bool);
/* increment */
bool limb_inc(limb_t *, limb_t);
/* decrement */
bool limb_dec(limb_t *, limb_t);
/* increment if carry, else copy */
bool limb_incc(limb_t *, limb_t, bool);
/* multiply and add, with two addends */
limb_t limb_fmaa(limb_t *, limb_t, limb_t, limb_t, limb_t);
/* divide, returning quotient and remainder */
limb_t limb_div(limb_t *, limb_t, limb_t, limb_t);
/* index of most significant nonzero digit */
Py_ssize_t limb_dsr(limb_t);
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
/* select high digits */
limb_t limb_high(limb_t, Py_ssize_t);
char limb_getdigit(limb_t, Py_ssize_t);
limb_t limb_setdigit(limb_t, Py_ssize_t, char);

bool limb_eq(limb_t, limb_t);
bool limb_le(limb_t, limb_t);
bool limb_lt(limb_t, limb_t);


limb_t limb_from_ulong(unsigned long *, unsigned long);
bool limb_to_ulong(unsigned long *, unsigned long, limb_t);
unsigned long limb_hash(limb_t);
