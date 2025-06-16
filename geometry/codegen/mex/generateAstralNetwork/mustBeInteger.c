/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mustBeInteger.c
 *
 * Code generation for function 'mustBeInteger'
 *
 */

/* Include files */
#include "mustBeInteger.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRTEInfo f_emlrtRTEI =
    {
        13,              /* lineNo */
        27,              /* colNo */
        "mustBeInteger", /* fName */
        "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/validators/"
        "mustBeInteger.m" /* pName */
};

/* Function Definitions */
void mustBeInteger(const emlrtStack *sp, const real_T A[2])
{
  int32_T k;
  boolean_T p;
  p = true;
  for (k = 0; k < 2; k++) {
    if (p) {
      real_T d;
      d = A[k];
      if (muDoubleScalarIsInf(d) || muDoubleScalarIsNaN(d) ||
          (!(d == muDoubleScalarFloor(d)))) {
        p = false;
      }
    } else {
      p = false;
    }
  }
  if (!p) {
    emlrtErrorWithMessageIdR2018a(sp, &f_emlrtRTEI,
                                  "MATLAB:validators:mustBeInteger",
                                  "MATLAB:validators:mustBeInteger", 0);
  }
}

/* End of code generation (mustBeInteger.c) */
