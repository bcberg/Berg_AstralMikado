/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "eml_int_forloop_overflow_check.h"
#include "generateAstralNetwork_data.h"
#include "generateAstralNetwork_emxutil.h"
#include "generateAstralNetwork_types.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include <emmintrin.h>

/* Variable Definitions */
static emlrtRSInfo bb_emlrtRSI = {
    20,    /* lineNo */
    "sum", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/datafun/sum.m" /* pathName
                                                                     */
};

static emlrtRSInfo cb_emlrtRSI =
    {
        99,        /* lineNo */
        "sumprod", /* fcnName */
        "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/datafun/private/"
        "sumprod.m" /* pathName */
};

static emlrtRSInfo db_emlrtRSI = {
    86,                      /* lineNo */
    "combineVectorElements", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/datafun/private/"
    "combineVectorElements.m" /* pathName */
};

static emlrtRSInfo eb_emlrtRSI = {
    112,                /* lineNo */
    "blockedSummation", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRSInfo fb_emlrtRSI = {
    173,                /* lineNo */
    "colMajorFlatIter", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRSInfo gb_emlrtRSI = {
    192,                /* lineNo */
    "colMajorFlatIter", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRTEInfo ec_emlrtRTEI = {
    146,                /* lineNo */
    24,                 /* colNo */
    "blockedSummation", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pName */
};

/* Function Definitions */
void sum(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y)
{
  jmp_buf *volatile emlrtJBStack;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  const real_T *x_data;
  real_T *y_data;
  int32_T sum_numThreads;
  int32_T xj;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  x_data = x->data;
  st.site = &bb_emlrtRSI;
  b_st.site = &cb_emlrtRSI;
  c_st.site = &db_emlrtRSI;
  if (x->size[0] == 0) {
    y->size[0] = 0;
  } else {
    int32_T scalarLB;
    int32_T vstride;
    int32_T vstride_tmp;
    d_st.site = &eb_emlrtRSI;
    vstride_tmp = x->size[0];
    vstride = x->size[0];
    scalarLB = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real_T(&d_st, y, scalarLB, &ec_emlrtRTEI);
    y_data = y->data;
    e_st.site = &fb_emlrtRSI;
    if (x->size[0] > 2147483646) {
      f_st.site = &t_emlrtRSI;
      check_forloop_overflow_error(&f_st);
    }
    if (x->size[0] < 3200) {
      for (xj = 0; xj < vstride_tmp; xj++) {
        y_data[xj] = x_data[xj];
      }
    } else {
      emlrtEnterParallelRegion(&d_st, omp_in_parallel());
      emlrtPushJmpBuf(&d_st, &emlrtJBStack);
      sum_numThreads =
          emlrtAllocRegionTLSs(d_st.tls, omp_in_parallel(),
                               omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(sum_numThreads)

      for (xj = 0; xj < vstride; xj++) {
        y_data[xj] = x_data[xj];
      }
      emlrtPopJmpBuf(&d_st, &emlrtJBStack);
      emlrtExitParallelRegion(&d_st, omp_in_parallel());
    }
    e_st.site = &gb_emlrtRSI;
    scalarLB = (vstride / 2) << 1;
    vstride_tmp = scalarLB - 2;
    for (xj = 0; xj <= vstride_tmp; xj += 2) {
      __m128d r;
      r = _mm_loadu_pd(&y_data[xj]);
      _mm_storeu_pd(&y_data[xj],
                    _mm_add_pd(r, _mm_loadu_pd(&x_data[vstride + xj])));
    }
    if (vstride - scalarLB < 3200) {
      for (xj = scalarLB; xj < vstride; xj++) {
        y_data[xj] += x_data[vstride + xj];
      }
    } else {
      emlrtEnterParallelRegion(&d_st, omp_in_parallel());
      emlrtPushJmpBuf(&d_st, &emlrtJBStack);
      sum_numThreads =
          emlrtAllocRegionTLSs(d_st.tls, omp_in_parallel(),
                               omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(sum_numThreads)

      for (xj = scalarLB; xj < vstride; xj++) {
        y_data[xj] += x_data[vstride + xj];
      }
      emlrtPopJmpBuf(&d_st, &emlrtJBStack);
      emlrtExitParallelRegion(&d_st, omp_in_parallel());
    }
  }
}

/* End of code generation (sum.c) */
