/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mod.c
 *
 * Code generation for function 'mod'
 *
 */

/* Include files */
#include "mod.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
real_T b_mod(real_T x, real_T y)
{
  real_T r;
  r = x;
  if (muDoubleScalarIsNaN(x) || muDoubleScalarIsNaN(y) ||
      muDoubleScalarIsInf(x)) {
    r = rtNaN;
  } else if (!muDoubleScalarIsInf(y)) {
    r = muDoubleScalarRem(x, y);
    if (r == 0.0) {
      r = 0.0;
    }
  }
  return r;
}

/* End of code generation (mod.c) */
