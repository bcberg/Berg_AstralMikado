/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rand.c
 *
 * Code generation for function 'rand'
 *
 */

/* Include files */
#include "rand.h"
#include "generateAstralNetwork_emxutil.h"
#include "generateAstralNetwork_types.h"
#include "mustBeInteger.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo g_emlrtRSI = {
    103,    /* lineNo */
    "rand", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/randfun/rand.m" /* pathName
                                                                      */
};

static emlrtRSInfo h_emlrtRSI = {
    104,    /* lineNo */
    "rand", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/randfun/rand.m" /* pathName
                                                                      */
};

static emlrtRSInfo i_emlrtRSI = {
    20,                /* lineNo */
    "zerosWithChecks", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/"
    "zerosWithChecks.m" /* pathName */
};

static emlrtDCInfo f_emlrtDCI = {
    23,                /* lineNo */
    13,                /* colNo */
    "zerosWithChecks", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/"
    "zerosWithChecks.m", /* pName */
    4                    /* checkKind */
};

static emlrtRTEInfo cb_emlrtRTEI = {
    103,    /* lineNo */
    24,     /* colNo */
    "rand", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/randfun/rand.m" /* pName */
};

/* Function Definitions */
void b_rand(const emlrtStack *sp, const real_T varargin_1[2],
            emxArray_real_T *r)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T *r_data;
  int32_T i;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &g_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  b_st.site = &i_emlrtRSI;
  mustBeInteger(&b_st, varargin_1);
  if (!(varargin_1[0] >= 0.0)) {
    emlrtNonNegativeCheckR2012b(varargin_1[0], &f_emlrtDCI, &st);
  }
  if (!(varargin_1[1] >= 0.0)) {
    emlrtNonNegativeCheckR2012b(varargin_1[1], &f_emlrtDCI, &st);
  }
  i = r->size[0] * r->size[1];
  r->size[0] = (int32_T)varargin_1[0];
  r->size[1] = 2;
  emxEnsureCapacity_real_T(&st, r, i, &cb_emlrtRTEI);
  r_data = r->data;
  st.site = &h_emlrtRSI;
  if ((int32_T)varargin_1[0] != 0) {
    emlrtRandu(&r_data[0], (int32_T)varargin_1[0] << 1);
  }
}

void c_rand(const emlrtStack *sp, const real_T varargin_1[2],
            emxArray_real_T *r)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T *r_data;
  int32_T i;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &g_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  b_st.site = &i_emlrtRSI;
  mustBeInteger(&b_st, varargin_1);
  if (!(varargin_1[0] >= 0.0)) {
    emlrtNonNegativeCheckR2012b(varargin_1[0], &f_emlrtDCI, &st);
  }
  if (!(varargin_1[1] >= 0.0)) {
    emlrtNonNegativeCheckR2012b(varargin_1[1], &f_emlrtDCI, &st);
  }
  i = r->size[0] * r->size[1];
  r->size[0] = (int32_T)varargin_1[0];
  r->size[1] = (int32_T)varargin_1[1];
  emxEnsureCapacity_real_T(&st, r, i, &cb_emlrtRTEI);
  r_data = r->data;
  st.site = &h_emlrtRSI;
  if (((int32_T)varargin_1[0] != 0) && ((int32_T)varargin_1[1] != 0)) {
    emlrtRandu(&r_data[0], (int32_T)varargin_1[0] * (int32_T)varargin_1[1]);
  }
}

/* End of code generation (rand.c) */
