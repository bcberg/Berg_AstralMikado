/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 */

/* Include files */
#include "repmat.h"
#include "assertValidSizeArg.h"
#include "eml_int_forloop_overflow_check.h"
#include "generateAstralNetwork_data.h"
#include "generateAstralNetwork_emxutil.h"
#include "generateAstralNetwork_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo y_emlrtRSI = {
    38,       /* lineNo */
    "repmat", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/repmat.m" /* pathName
                                                                      */
};

static emlrtRSInfo ab_emlrtRSI = {
    79,       /* lineNo */
    "repmat", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/repmat.m" /* pathName
                                                                      */
};

static emlrtRTEInfo dc_emlrtRTEI = {
    69,       /* lineNo */
    28,       /* colNo */
    "repmat", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/repmat.m" /* pName */
};

/* Function Definitions */
void repmat(const emlrtStack *sp, const real_T a[2], const real_T varargin_1[2],
            emxArray_real_T *b)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T *b_data;
  int32_T i;
  int32_T ibmat;
  int32_T itilerow;
  int32_T jcol;
  boolean_T overflow;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &y_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  assertValidSizeArg(&st, varargin_1);
  i = (int32_T)varargin_1[0];
  ibmat = b->size[0] * b->size[1];
  b->size[0] = (int32_T)varargin_1[0];
  b->size[1] = 2;
  emxEnsureCapacity_real_T(sp, b, ibmat, &dc_emlrtRTEI);
  b_data = b->data;
  overflow = ((int32_T)varargin_1[0] > 2147483646);
  for (jcol = 0; jcol < 2; jcol++) {
    ibmat = jcol * i;
    st.site = &ab_emlrtRSI;
    if (overflow) {
      b_st.site = &t_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }
    for (itilerow = 0; itilerow < i; itilerow++) {
      b_data[ibmat + itilerow] = a[jcol];
    }
  }
}

/* End of code generation (repmat.c) */
