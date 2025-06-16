/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sort.c
 *
 * Code generation for function 'sort'
 *
 */

/* Include files */
#include "sort.h"
#include "eml_int_forloop_overflow_check.h"
#include "generateAstralNetwork_data.h"
#include "generateAstralNetwork_emxutil.h"
#include "generateAstralNetwork_types.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"

/* Variable Definitions */
static emlrtRSInfo ic_emlrtRSI = {
    76,     /* lineNo */
    "sort", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pathName
                                                                        */
};

static emlrtRSInfo jc_emlrtRSI = {
    79,     /* lineNo */
    "sort", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pathName
                                                                        */
};

static emlrtRSInfo kc_emlrtRSI = {
    81,     /* lineNo */
    "sort", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pathName
                                                                        */
};

static emlrtRSInfo lc_emlrtRSI = {
    84,     /* lineNo */
    "sort", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pathName
                                                                        */
};

static emlrtRSInfo mc_emlrtRSI = {
    87,     /* lineNo */
    "sort", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pathName
                                                                        */
};

static emlrtRSInfo nc_emlrtRSI = {
    90,     /* lineNo */
    "sort", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pathName
                                                                        */
};

static emlrtRTEInfo uc_emlrtRTEI = {
    56,     /* lineNo */
    24,     /* colNo */
    "sort", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pName
                                                                        */
};

static emlrtRTEInfo vc_emlrtRTEI = {
    75,     /* lineNo */
    26,     /* colNo */
    "sort", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pName
                                                                        */
};

static emlrtRTEInfo wc_emlrtRTEI = {
    56,     /* lineNo */
    1,      /* colNo */
    "sort", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pName
                                                                        */
};

static emlrtRTEInfo xc_emlrtRTEI = {
    1,      /* lineNo */
    20,     /* colNo */
    "sort", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/sort.m" /* pName
                                                                        */
};

/* Function Definitions */
void sort(const emlrtStack *sp, emxArray_real_T *x, emxArray_int32_T *idx)
{
  emlrtStack b_st;
  emlrtStack st;
  emxArray_int32_T *iidx;
  emxArray_real_T *vwork;
  real_T *vwork_data;
  real_T *x_data;
  int32_T i;
  int32_T k;
  int32_T vlen;
  int32_T *idx_data;
  int32_T *iidx_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  x_data = x->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  k = x->size[0];
  vlen = x->size[0] - 1;
  emxInit_real_T(sp, &vwork, 1, &wc_emlrtRTEI);
  i = vwork->size[0];
  vwork->size[0] = k;
  emxEnsureCapacity_real_T(sp, vwork, i, &uc_emlrtRTEI);
  vwork_data = vwork->data;
  i = idx->size[0];
  idx->size[0] = k;
  emxEnsureCapacity_int32_T(sp, idx, i, &vc_emlrtRTEI);
  idx_data = idx->data;
  st.site = &ic_emlrtRSI;
  st.site = &jc_emlrtRSI;
  st.site = &kc_emlrtRSI;
  st.site = &lc_emlrtRSI;
  if (k > 2147483646) {
    b_st.site = &t_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }
  for (k = 0; k <= vlen; k++) {
    vwork_data[k] = x_data[k];
  }
  emxInit_int32_T(sp, &iidx, &xc_emlrtRTEI);
  st.site = &mc_emlrtRSI;
  sortIdx(&st, vwork, iidx);
  iidx_data = iidx->data;
  vwork_data = vwork->data;
  st.site = &nc_emlrtRSI;
  for (k = 0; k <= vlen; k++) {
    x_data[k] = vwork_data[k];
    idx_data[k] = iidx_data[k];
  }
  emxFree_int32_T(sp, &iidx);
  emxFree_real_T(sp, &vwork);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (sort.c) */
