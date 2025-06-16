/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * error.c
 *
 * Code generation for function 'error'
 *
 */

/* Include files */
#include "error.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtMCInfo emlrtMCI = {
    27,                                                            /* lineNo */
    5,                                                             /* colNo */
    "error",                                                       /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/lang/error.m" /* pName */
};

static emlrtRSInfo kd_emlrtRSI = {
    27,                                                            /* lineNo */
    "error",                                                       /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/lang/error.m" /* pathName
                                                                    */
};

/* Function Declarations */
static void c_error(const emlrtStack *sp, const mxArray *m,
                    emlrtMCInfo *location);

/* Function Definitions */
static void c_error(const emlrtStack *sp, const mxArray *m,
                    emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = m;
  emlrtCallMATLABR2012b((emlrtConstCTX)sp, 0, NULL, 1, &pArray, "error", true,
                        location);
}

void b_error(const emlrtStack *sp)
{
  static const int32_T iv[2] = {1, 89};
  static const char_T varargin_1[89] = {
      'N', 'e', 't', 'w', 'o', 'r', 'k', ' ', 'g', 'e', 'n', 'e', 'r', 'a', 't',
      'i', 'o', 'n', ' ', 'f', 'a', 'i', 'l', 'e', 'd', ' ', '(', 'n', 'u', 'm',
      'A', 's', 't', 'e', 'r', 's', ' ', '=', ' ', '0', ')', '.', ' ', 'A', 'r',
      'g', 'u', 'm', 'e', 'n', 't', ' ', 'n', 'u', 'm', 'A', 's', 't', 'e', 'r',
      's', ' ', 'm', 'u', 's', 't', ' ', 'b', 'e', ' ', 'a', ' ', 'p', 'o', 's',
      'i', 't', 'i', 'v', 'e', ' ', 'i', 'n', 't', 'e', 'g', 'e', 'r', '.'};
  emlrtStack st;
  const mxArray *m;
  const mxArray *y;
  st.prev = sp;
  st.tls = sp->tls;
  y = NULL;
  m = emlrtCreateCharArray(2, &iv[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 89, m, &varargin_1[0]);
  emlrtAssign(&y, m);
  st.site = &kd_emlrtRSI;
  c_error(&st, y, &emlrtMCI);
}

/* End of code generation (error.c) */
