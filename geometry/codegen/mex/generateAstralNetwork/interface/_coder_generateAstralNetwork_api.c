/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_generateAstralNetwork_api.c
 *
 * Code generation for function '_coder_generateAstralNetwork_api'
 *
 */

/* Include files */
#include "_coder_generateAstralNetwork_api.h"
#include "generateAstralNetwork.h"
#include "generateAstralNetwork_data.h"
#include "generateAstralNetwork_emxutil.h"
#include "generateAstralNetwork_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo tc_emlrtRTEI = {
    1,                                  /* lineNo */
    1,                                  /* colNo */
    "_coder_generateAstralNetwork_api", /* fName */
    ""                                  /* pName */
};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u);

static boolean_T c_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static const mxArray *c_emlrt_marshallOut(const struct1_T u);

static boolean_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static const mxArray *d_emlrt_marshallOut(const emxArray_real_T *u);

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static const mxArray *e_emlrt_marshallOut(const struct2_T u);

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier);

static const mxArray *emlrt_marshallOut(const struct0_T u);

static boolean_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *m;
  const mxArray *y;
  const real_T *u_data;
  real_T *pData;
  int32_T iv[2];
  int32_T b_i;
  int32_T i;
  u_data = u->data;
  y = NULL;
  iv[0] = u->size[0];
  iv[1] = 2;
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < u->size[0]; b_i++) {
    pData[i] = u_data[b_i];
    i++;
  }
  for (b_i = 0; b_i < u->size[0]; b_i++) {
    pData[i] = u_data[b_i + u->size[0]];
    i++;
  }
  emlrtAssign(&y, m);
  return y;
}

static boolean_T c_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  boolean_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *c_emlrt_marshallOut(const struct1_T u)
{
  static const char_T *sv[2] = {"centerCross", "filCross"};
  const mxArray *y;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 2, (const char_T **)&sv[0]));
  emlrtSetFieldR2017b(y, 0, "centerCross", d_emlrt_marshallOut(u.centerCross),
                      0);
  emlrtSetFieldR2017b(y, 0, "filCross", b_emlrt_marshallOut(u.filCross), 1);
  return y;
}

static boolean_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId)
{
  boolean_T y;
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *d_emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *m;
  const mxArray *y;
  const real_T *u_data;
  real_T *pData;
  int32_T iv[2];
  int32_T b_i;
  int32_T c_i;
  int32_T i;
  u_data = u->data;
  y = NULL;
  iv[0] = u->size[0];
  iv[1] = u->size[1];
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < u->size[1]; b_i++) {
    for (c_i = 0; c_i < u->size[0]; c_i++) {
      pData[i] = u_data[c_i + u->size[0] * b_i];
      i++;
    }
  }
  emlrtAssign(&y, m);
  return y;
}

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static const mxArray *e_emlrt_marshallOut(const struct2_T u)
{
  static const char_T *sv[2] = {"centers", "orients"};
  const mxArray *y;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 2, (const char_T **)&sv[0]));
  emlrtSetFieldR2017b(y, 0, "centers", b_emlrt_marshallOut(u.centers), 0);
  emlrtSetFieldR2017b(y, 0, "orients", d_emlrt_marshallOut(u.orients), 1);
  return y;
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *emlrt_marshallOut(const struct0_T u)
{
  static const char_T *sv[4] = {"nodes", "augNodes", "springs", "ends"};
  const mxArray *b_y;
  const mxArray *m;
  const mxArray *y;
  real_T *pData;
  int32_T iv[2];
  int32_T b_i;
  int32_T i;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 4, (const char_T **)&sv[0]));
  emlrtSetFieldR2017b(y, 0, "nodes", b_emlrt_marshallOut(u.nodes), 0);
  emlrtSetFieldR2017b(y, 0, "augNodes", b_emlrt_marshallOut(u.augNodes), 1);
  b_y = NULL;
  iv[0] = u.springs->size[0];
  iv[1] = 4;
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < u.springs->size[0]; b_i++) {
    pData[i] = u.springs->data[b_i];
    i++;
  }
  for (b_i = 0; b_i < u.springs->size[0]; b_i++) {
    pData[i] = u.springs->data[b_i + u.springs->size[0]];
    i++;
  }
  for (b_i = 0; b_i < u.springs->size[0]; b_i++) {
    pData[i] = u.springs->data[b_i + u.springs->size[0] * 2];
    i++;
  }
  for (b_i = 0; b_i < u.springs->size[0]; b_i++) {
    pData[i] = u.springs->data[b_i + u.springs->size[0] * 3];
    i++;
  }
  emlrtAssign(&b_y, m);
  emlrtSetFieldR2017b(y, 0, "springs", b_y, 2);
  emlrtSetFieldR2017b(y, 0, "ends", b_emlrt_marshallOut(u.ends), 3);
  return y;
}

static boolean_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  boolean_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "logical", false, 0U,
                          (const void *)&dims);
  ret = *emlrtMxGetLogicals(src);
  emlrtDestroyArray(&src);
  return ret;
}

void generateAstralNetwork_api(const mxArray *const prhs[5], int32_T nlhs,
                               const mxArray *plhs[3])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  struct0_T network;
  struct1_T crossings;
  struct2_T asters;
  real_T D;
  real_T astralNum;
  real_T l;
  real_T numAsters;
  boolean_T nodesOnly;
  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  /* Marshall function inputs */
  numAsters = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "numAsters");
  l = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "l");
  D = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "D");
  astralNum = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "astralNum");
  nodesOnly = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "nodesOnly");
  /* Invoke the target function */
  emxInitStruct_struct0_T(&st, &network, &tc_emlrtRTEI);
  emxInitStruct_struct1_T(&st, &crossings, &tc_emlrtRTEI);
  emxInitStruct_struct2_T(&st, &asters, &tc_emlrtRTEI);
  generateAstralNetwork(&st, numAsters, l, D, astralNum, nodesOnly, &network,
                        &crossings, &asters);
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(network);
  emxFreeStruct_struct0_T(&st, &network);
  if (nlhs > 1) {
    plhs[1] = c_emlrt_marshallOut(crossings);
  }
  emxFreeStruct_struct1_T(&st, &crossings);
  if (nlhs > 2) {
    plhs[2] = e_emlrt_marshallOut(asters);
  }
  emxFreeStruct_struct2_T(&st, &asters);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_generateAstralNetwork_api.c) */
