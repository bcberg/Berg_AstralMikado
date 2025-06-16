/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * colon.c
 *
 * Code generation for function 'colon'
 *
 */

/* Include files */
#include "colon.h"
#include "generateAstralNetwork_emxutil.h"
#include "generateAstralNetwork_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRTEInfo cc_emlrtRTEI = {
    320,                                                          /* lineNo */
    20,                                                           /* colNo */
    "colon",                                                      /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/ops/colon.m" /* pName */
};

/* Function Definitions */
void eml_float_colon(const emlrtStack *sp, real_T a, real_T b,
                     emxArray_real_T *y)
{
  real_T apnd;
  real_T cdiff;
  real_T ndbl;
  real_T *y_data;
  int32_T k;
  int32_T n;
  int32_T nm1d2;
  ndbl = muDoubleScalarFloor((b - a) + 0.5);
  apnd = a + ndbl;
  cdiff = apnd - b;
  if (muDoubleScalarAbs(cdiff) <
      4.4408920985006262E-16 *
          muDoubleScalarMax(muDoubleScalarAbs(a), muDoubleScalarAbs(b))) {
    ndbl++;
    apnd = b;
  } else if (cdiff > 0.0) {
    apnd = a + (ndbl - 1.0);
  } else {
    ndbl++;
  }
  if (ndbl >= 0.0) {
    n = (int32_T)ndbl;
  } else {
    n = 0;
  }
  nm1d2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity_real_T(sp, y, nm1d2, &cc_emlrtRTEI);
  y_data = y->data;
  if (n > 0) {
    y_data[0] = a;
    if (n > 1) {
      y_data[n - 1] = apnd;
      nm1d2 = (n - 1) / 2;
      for (k = 0; k <= nm1d2 - 2; k++) {
        y_data[k + 1] = a + ((real_T)k + 1.0);
        y_data[(n - k) - 2] = apnd - ((real_T)k + 1.0);
      }
      if (nm1d2 << 1 == n - 1) {
        y_data[nm1d2] = (a + apnd) / 2.0;
      } else {
        y_data[nm1d2] = a + (real_T)nm1d2;
        y_data[nm1d2 + 1] = apnd - (real_T)nm1d2;
      }
    }
  }
}

/* End of code generation (colon.c) */
