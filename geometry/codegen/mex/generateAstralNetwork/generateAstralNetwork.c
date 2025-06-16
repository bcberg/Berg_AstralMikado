/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * generateAstralNetwork.c
 *
 * Code generation for function 'generateAstralNetwork'
 *
 */

/* Include files */
#include "generateAstralNetwork.h"
#include "assertValidSizeArg.h"
#include "colon.h"
#include "eml_int_forloop_overflow_check.h"
#include "error.h"
#include "generateAstralNetwork_data.h"
#include "generateAstralNetwork_emxutil.h"
#include "generateAstralNetwork_types.h"
#include "indexShapeCheck.h"
#include "mod.h"
#include "rand.h"
#include "repmat.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "sum.h"
#include "sumMatrixIncludeNaN.h"
#include "warning.h"
#include "mwmathutil.h"
#include "omp.h"
#include <emmintrin.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    30,                                                            /* lineNo */
    "generateAstralNetwork",                                       /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo b_emlrtRSI = {
    33,                                                            /* lineNo */
    "generateAstralNetwork",                                       /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo c_emlrtRSI = {
    34,                                                            /* lineNo */
    "generateAstralNetwork",                                       /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo d_emlrtRSI = {
    35,                                                            /* lineNo */
    "generateAstralNetwork",                                       /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo e_emlrtRSI = {
    38,                                                            /* lineNo */
    "generateAstralNetwork",                                       /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo f_emlrtRSI = {
    48,                                                            /* lineNo */
    "generateAstralNetwork",                                       /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo j_emlrtRSI = {
    136,                                                           /* lineNo */
    "findNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo k_emlrtRSI = {
    134,                                                           /* lineNo */
    "findNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo l_emlrtRSI = {
    133,                                                           /* lineNo */
    "findNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo m_emlrtRSI = {
    93,                                                            /* lineNo */
    "findNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo n_emlrtRSI = {
    91,                                                            /* lineNo */
    "findNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo o_emlrtRSI = {
    90,                                                            /* lineNo */
    "findNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo p_emlrtRSI = {
    81,                                                            /* lineNo */
    "findNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo q_emlrtRSI = {
    80,                                                            /* lineNo */
    "findNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo r_emlrtRSI = {
    11,                                                           /* lineNo */
    "cos",                                                        /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elfun/cos.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI = {
    33,                           /* lineNo */
    "applyScalarFunctionInPlace", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/"
    "applyScalarFunctionInPlace.m" /* pathName */
};

static emlrtRSInfo u_emlrtRSI = {
    11,                                                           /* lineNo */
    "sin",                                                        /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elfun/sin.m" /* pathName */
};

static emlrtRSInfo v_emlrtRSI = {
    28,                                                           /* lineNo */
    "colon",                                                      /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/ops/colon.m" /* pathName */
};

static emlrtRSInfo w_emlrtRSI = {
    125,                                                          /* lineNo */
    "colon",                                                      /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/ops/colon.m" /* pathName */
};

static emlrtRSInfo hb_emlrtRSI = {
    44,       /* lineNo */
    "mpower", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/matfun/mpower.m" /* pathName
                                                                       */
};

static emlrtRSInfo ib_emlrtRSI = {
    71,                                                           /* lineNo */
    "power",                                                      /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/ops/power.m" /* pathName */
};

static emlrtRSInfo jb_emlrtRSI = {
    40,                  /* lineNo */
    "reshapeSizeChecks", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtRSInfo kb_emlrtRSI = {
    343,                                                           /* lineNo */
    "defineSprings",                                               /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo lb_emlrtRSI = {
    334,                                                           /* lineNo */
    "defineSprings",                                               /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo nb_emlrtRSI = {
    318,                                                           /* lineNo */
    "defineSprings",                                               /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo ob_emlrtRSI = {
    308,                                                           /* lineNo */
    "defineSprings",                                               /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo pb_emlrtRSI = {
    300,                                                           /* lineNo */
    "defineSprings",                                               /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo qb_emlrtRSI = {
    293,                                                           /* lineNo */
    "defineSprings",                                               /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo rb_emlrtRSI = {
    280,                                                           /* lineNo */
    "defineSprings",                                               /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo sb_emlrtRSI = {
    271,                                                           /* lineNo */
    "defineSprings",                                               /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo tb_emlrtRSI = {
    207,                                                           /* lineNo */
    "sortNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo ub_emlrtRSI = {
    220,                                                           /* lineNo */
    "sortNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo vb_emlrtRSI = {
    222,                                                           /* lineNo */
    "sortNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo wb_emlrtRSI = {
    223,                                                           /* lineNo */
    "sortNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo xb_emlrtRSI = {
    226,                                                           /* lineNo */
    "sortNodes",                                                   /* fcnName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pathName
                                                                    */
};

static emlrtRSInfo yb_emlrtRSI = {
    35,                                                            /* lineNo */
    "find",                                                        /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pathName
                                                                    */
};

static emlrtRSInfo ac_emlrtRSI = {
    158,                                                           /* lineNo */
    "eml_find",                                                    /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pathName
                                                                    */
};

static emlrtRSInfo bc_emlrtRSI = {
    138,                                                           /* lineNo */
    "eml_find",                                                    /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pathName
                                                                    */
};

static emlrtRSInfo cc_emlrtRSI = {
    253,                                                           /* lineNo */
    "find_first_nonempty_triples",                                 /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pathName
                                                                    */
};

static emlrtRSInfo dc_emlrtRSI = {
    254,                                                           /* lineNo */
    "find_first_nonempty_triples",                                 /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pathName
                                                                    */
};

static emlrtRSInfo ec_emlrtRSI = {
    255,                                                           /* lineNo */
    "find_first_nonempty_triples",                                 /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pathName
                                                                    */
};

static emlrtRSInfo gc_emlrtRSI = {
    70,    /* lineNo */
    "dot", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/specfun/dot.m" /* pathName
                                                                     */
};

static emlrtRSInfo hc_emlrtRSI = {
    37,     /* lineNo */
    "sort", /* fcnName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/datafun/sort.m" /* pathName
                                                                      */
};

static emlrtBCInfo emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    321,                                                            /* lineNo */
    18,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    347,                                                            /* lineNo */
    18,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    301,                                                            /* lineNo */
    18,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    301,                                                            /* lineNo */
    54,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    300,                                                            /* lineNo */
    18,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    281,                                                            /* lineNo */
    18,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    281,                                                            /* lineNo */
    36,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    280,                                                            /* lineNo */
    18,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    275,                                                            /* lineNo */
    18,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    312,                                                            /* lineNo */
    18,                                                             /* colNo */
    "ends",                                                         /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtDCInfo emlrtDCI = {
    264,                                                            /* lineNo */
    1,                                                              /* colNo */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = {
    263,                                                            /* lineNo */
    1,                                                              /* colNo */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = {
    264,                                                            /* lineNo */
    14,                                                             /* colNo */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    4 /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = {
    264,                                                            /* lineNo */
    14,                                                             /* colNo */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = {
    263,                                                            /* lineNo */
    17,                                                             /* colNo */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtRTEInfo emlrtRTEI = {
    74,                  /* lineNo */
    13,                  /* colNo */
    "reshapeSizeChecks", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pName */
};

static emlrtRTEInfo b_emlrtRTEI = {
    79,                  /* lineNo */
    23,                  /* colNo */
    "reshapeSizeChecks", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pName */
};

static emlrtRTEInfo c_emlrtRTEI = {
    81,                  /* lineNo */
    23,                  /* colNo */
    "reshapeSizeChecks", /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pName */
};

static emlrtBCInfo k_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    354,                                                            /* lineNo */
    25,                                                             /* colNo */
    "springs",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    354,                                                            /* lineNo */
    23,                                                             /* colNo */
    "springs",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtECInfo emlrtECI = {
    -1,                                                            /* nDims */
    320,                                                           /* lineNo */
    13,                                                            /* colNo */
    "defineSprings",                                               /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo m_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    320,                                                            /* lineNo */
    21,                                                             /* colNo */
    "springs",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    342,                                                            /* lineNo */
    30,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo o_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    333,                                                            /* lineNo */
    31,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo p_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    332,                                                            /* lineNo */
    38,                                                             /* colNo */
    "thisOrder",                                                    /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo q_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    331,                                                            /* lineNo */
    31,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo r_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    330,                                                            /* lineNo */
    38,                                                             /* colNo */
    "thisOrder",                                                    /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo s_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    324,                                                            /* lineNo */
    26,                                                             /* colNo */
    "centers",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo t_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    317,                                                            /* lineNo */
    32,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo u_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    316,                                                            /* lineNo */
    26,                                                             /* colNo */
    "centers",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtRTEInfo d_emlrtRTEI = {
    307,                                                           /* lineNo */
    15,                                                            /* colNo */
    "defineSprings",                                               /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo v_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    299,                                                            /* lineNo */
    30,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo w_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    292,                                                            /* lineNo */
    31,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo x_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    291,                                                            /* lineNo */
    38,                                                             /* colNo */
    "thisOrder",                                                    /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo y_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    290,                                                            /* lineNo */
    31,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo ab_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    289,                                                            /* lineNo */
    38,                                                             /* colNo */
    "thisOrder",                                                    /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo bb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    284,                                                            /* lineNo */
    26,                                                             /* colNo */
    "centers",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo cb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    279,                                                            /* lineNo */
    32,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo db_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    278,                                                            /* lineNo */
    26,                                                             /* colNo */
    "centers",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtRTEInfo e_emlrtRTEI = {
    270,                                                           /* lineNo */
    15,                                                            /* colNo */
    "defineSprings",                                               /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo eb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    345,                                                            /* lineNo */
    21,                                                             /* colNo */
    "springs",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo fb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    296,                                                            /* lineNo */
    25,                                                             /* colNo */
    "springs",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo gb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    337,                                                            /* lineNo */
    25,                                                             /* colNo */
    "springs",                                                      /* aName */
    "defineSprings",                                                /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo hb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    188,                                                            /* lineNo */
    27,                                                             /* colNo */
    "filCross",                                                     /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo ib_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    188,                                                            /* lineNo */
    25,                                                             /* colNo */
    "filCross",                                                     /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo jb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    187,                                                            /* lineNo */
    21,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo kb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    187,                                                            /* lineNo */
    19,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo lb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    178,                                                            /* lineNo */
    26,                                                             /* colNo */
    "filCross",                                                     /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo mb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    177,                                                            /* lineNo */
    23,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo nb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    176,                                                            /* lineNo */
    66,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo ob_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    176,                                                            /* lineNo */
    57,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo pb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    176,                                                            /* lineNo */
    33,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo qb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    175,                                                            /* lineNo */
    68,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo rb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    175,                                                            /* lineNo */
    59,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtECInfo b_emlrtECI = {
    -1,                                                            /* nDims */
    173,                                                           /* lineNo */
    21,                                                            /* colNo */
    "findNodes",                                                   /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo sb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    173,                                                            /* lineNo */
    33,                                                             /* colNo */
    "filCross",                                                     /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo tb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    173,                                                            /* lineNo */
    30,                                                             /* colNo */
    "filCross",                                                     /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtECInfo c_emlrtECI = {
    -1,                                                            /* nDims */
    170,                                                           /* lineNo */
    21,                                                            /* colNo */
    "findNodes",                                                   /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo ub_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    170,                                                            /* lineNo */
    30,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo vb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    170,                                                            /* lineNo */
    27,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo wb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    163,                                                            /* lineNo */
    32,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo xb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    163,                                                            /* lineNo */
    23,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo yb_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    162,                                                            /* lineNo */
    36,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtRTEInfo g_emlrtRTEI = {
    138,                                                           /* lineNo */
    19,                                                            /* colNo */
    "findNodes",                                                   /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtECInfo d_emlrtECI = {
    1,                                                             /* nDims */
    134,                                                           /* lineNo */
    29,                                                            /* colNo */
    "findNodes",                                                   /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo ac_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    134,                                                            /* lineNo */
    44,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtRTEInfo h_emlrtRTEI = {
    131,                                                           /* lineNo */
    15,                                                            /* colNo */
    "findNodes",                                                   /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo bc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    124,                                                            /* lineNo */
    26,                                                             /* colNo */
    "filCross",                                                     /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo cc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    123,                                                            /* lineNo */
    23,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo dc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    122,                                                            /* lineNo */
    52,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo ec_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    122,                                                            /* lineNo */
    33,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo fc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    121,                                                            /* lineNo */
    54,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo gc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    121,                                                            /* lineNo */
    33,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtECInfo e_emlrtECI = {
    -1,                                                            /* nDims */
    119,                                                           /* lineNo */
    21,                                                            /* colNo */
    "findNodes",                                                   /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo hc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    119,                                                            /* lineNo */
    33,                                                             /* colNo */
    "filCross",                                                     /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo ic_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    119,                                                            /* lineNo */
    30,                                                             /* colNo */
    "filCross",                                                     /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtECInfo f_emlrtECI = {
    -1,                                                            /* nDims */
    116,                                                           /* lineNo */
    21,                                                            /* colNo */
    "findNodes",                                                   /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo jc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    116,                                                            /* lineNo */
    30,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo kc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    116,                                                            /* lineNo */
    27,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo lc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    109,                                                            /* lineNo */
    27,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo mc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    95,                                                             /* lineNo */
    36,                                                             /* colNo */
    "closeEnough",                                                  /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtECInfo g_emlrtECI = {
    1,                                                             /* nDims */
    91,                                                            /* lineNo */
    29,                                                            /* colNo */
    "findNodes",                                                   /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtBCInfo nc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    91,                                                             /* lineNo */
    44,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtDCInfo g_emlrtDCI = {
    83,                                                             /* lineNo */
    15,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo h_emlrtDCI = {
    83,                                                             /* lineNo */
    15,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    4 /* checkKind */
};

static emlrtDCInfo i_emlrtDCI = {
    84,                                                             /* lineNo */
    18,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo j_emlrtDCI = {
    115,                                                            /* lineNo */
    35,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo k_emlrtDCI = {
    118,                                                            /* lineNo */
    38,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo l_emlrtDCI = {
    169,                                                            /* lineNo */
    35,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo m_emlrtDCI = {
    172,                                                            /* lineNo */
    38,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo n_emlrtDCI = {
    83,                                                             /* lineNo */
    1,                                                              /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo o_emlrtDCI = {
    84,                                                             /* lineNo */
    1,                                                              /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtBCInfo oc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    92,                                                             /* lineNo */
    21,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtDCInfo p_emlrtDCI = {
    92,                                                             /* lineNo */
    21,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtBCInfo pc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    135,                                                            /* lineNo */
    47,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtDCInfo q_emlrtDCI = {
    135,                                                            /* lineNo */
    47,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtBCInfo qc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    93,                                                             /* lineNo */
    33,                                                             /* colNo */
    "otherFils",                                                    /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo rc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    142,                                                            /* lineNo */
    23,                                                             /* colNo */
    "tooFar",                                                       /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo sc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    104,                                                            /* lineNo */
    29,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo tc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    104,                                                            /* lineNo */
    44,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo uc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    105,                                                            /* lineNo */
    23,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo vc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    105,                                                            /* lineNo */
    36,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo wc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    106,                                                            /* lineNo */
    27,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo xc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    106,                                                            /* lineNo */
    49,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo yc_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    107,                                                            /* lineNo */
    49,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo ad_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    108,                                                            /* lineNo */
    30,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo bd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    156,                                                            /* lineNo */
    29,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo cd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    156,                                                            /* lineNo */
    38,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtDCInfo r_emlrtDCI = {
    156,                                                            /* lineNo */
    38,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtBCInfo dd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    156,                                                            /* lineNo */
    57,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo ed_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    156,                                                            /* lineNo */
    66,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtDCInfo s_emlrtDCI = {
    156,                                                            /* lineNo */
    66,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtBCInfo fd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    157,                                                            /* lineNo */
    23,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo gd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    157,                                                            /* lineNo */
    32,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo hd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    157,                                                            /* lineNo */
    51,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo id_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    157,                                                            /* lineNo */
    60,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo jd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    158,                                                            /* lineNo */
    27,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo kd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    158,                                                            /* lineNo */
    49,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo ld_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    159,                                                            /* lineNo */
    49,                                                             /* colNo */
    "centers",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo md_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    160,                                                            /* lineNo */
    30,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo nd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    160,                                                            /* lineNo */
    39,                                                             /* colNo */
    "sines",                                                        /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo od_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    161,                                                            /* lineNo */
    28,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo pd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    161,                                                            /* lineNo */
    37,                                                             /* colNo */
    "cosines",                                                      /* aName */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtDCInfo t_emlrtDCI = {
    115,                                                            /* lineNo */
    21,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo u_emlrtDCI = {
    169,                                                            /* lineNo */
    21,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo v_emlrtDCI = {
    118,                                                            /* lineNo */
    21,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtDCInfo w_emlrtDCI = {
    172,                                                            /* lineNo */
    21,                                                             /* colNo */
    "findNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    1 /* checkKind */
};

static emlrtBCInfo qd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    218,                                                            /* lineNo */
    14,                                                             /* colNo */
    "centers",                                                      /* aName */
    "sortNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtECInfo h_emlrtECI = {
    1,                                                             /* nDims */
    222,                                                           /* lineNo */
    16,                                                            /* colNo */
    "sortNodes",                                                   /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo k_emlrtRTEI = {
    239,                                                           /* lineNo */
    1,                                                             /* colNo */
    "find_first_nonempty_triples",                                 /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pName */
};

static emlrtRTEInfo l_emlrtRTEI = {
    13,                                                            /* lineNo */
    9,                                                             /* colNo */
    "sqrt",                                                        /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elfun/sqrt.m" /* pName */
};

static emlrtRTEInfo m_emlrtRTEI = {
    12,                                                             /* lineNo */
    23,                                                             /* colNo */
    "dot",                                                          /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/specfun/dot.m" /* pName */
};

static emlrtBCInfo rd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    216,                                                            /* lineNo */
    25,                                                             /* colNo */
    "nodes",                                                        /* aName */
    "sortNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtBCInfo sd_emlrtBCI = {
    -1,                                                             /* iFirst */
    -1,                                                             /* iLast */
    227,                                                            /* lineNo */
    27,                                                             /* colNo */
    "nodesOnFil",                                                   /* aName */
    "sortNodes",                                                    /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m", /* pName */
    0 /* checkKind */
};

static emlrtRTEInfo p_emlrtRTEI = {
    33,                                                            /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo q_emlrtRTEI = {
    34,                                                            /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo r_emlrtRTEI = {
    38,                                                            /* lineNo */
    43,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo s_emlrtRTEI = {
    38,                                                            /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo t_emlrtRTEI = {
    263,                                                           /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo u_emlrtRTEI = {
    264,                                                           /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo v_emlrtRTEI = {
    51,                                                            /* lineNo */
    5,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo w_emlrtRTEI = {
    52,                                                            /* lineNo */
    5,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo x_emlrtRTEI = {
    1,                                                             /* lineNo */
    39,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo y_emlrtRTEI = {
    38,                                                            /* lineNo */
    35,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = {
    34,                                                            /* lineNo */
    27,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = {
    271,                                                           /* lineNo */
    9,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo db_emlrtRTEI = {
    80,                                                            /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = {
    81,                                                            /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo fb_emlrtRTEI = {
    83,                                                            /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo gb_emlrtRTEI = {
    84,                                                            /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo hb_emlrtRTEI = {
    187,                                                           /* lineNo */
    5,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo ib_emlrtRTEI = {
    90,                                                            /* lineNo */
    9,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo jb_emlrtRTEI = {
    133,                                                           /* lineNo */
    9,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo kb_emlrtRTEI = {
    188,                                                           /* lineNo */
    5,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo lb_emlrtRTEI = {
    92,                                                            /* lineNo */
    13,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo mb_emlrtRTEI = {
    135,                                                           /* lineNo */
    39,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo nb_emlrtRTEI =
    {
        31,            /* lineNo */
        30,            /* colNo */
        "unsafeSxfun", /* fName */
        "/usr/local/MATLAB/R2024a/toolbox/eml/eml/+coder/+internal/"
        "unsafeSxfun.m" /* pName */
};

static emlrtRTEInfo ob_emlrtRTEI = {
    59,                                                            /* lineNo */
    30,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = {
    136,                                                           /* lineNo */
    9,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo qb_emlrtRTEI = {
    114,                                                           /* lineNo */
    21,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo rb_emlrtRTEI = {
    115,                                                           /* lineNo */
    21,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo sb_emlrtRTEI = {
    117,                                                           /* lineNo */
    21,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo tb_emlrtRTEI = {
    118,                                                           /* lineNo */
    21,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo ub_emlrtRTEI = {
    168,                                                           /* lineNo */
    21,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo vb_emlrtRTEI = {
    169,                                                           /* lineNo */
    21,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo wb_emlrtRTEI = {
    171,                                                           /* lineNo */
    21,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo xb_emlrtRTEI = {
    172,                                                           /* lineNo */
    21,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo yb_emlrtRTEI = {
    91,                                                            /* lineNo */
    9,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo ac_emlrtRTEI = {
    134,                                                           /* lineNo */
    29,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo bc_emlrtRTEI = {
    93,                                                            /* lineNo */
    33,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo fc_emlrtRTEI = {
    207,                                                           /* lineNo */
    23,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo gc_emlrtRTEI = {
    195,                                                           /* lineNo */
    24,                                                            /* colNo */
    "find",                                                        /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pName */
};

static emlrtRTEInfo hc_emlrtRTEI = {
    36,                                                            /* lineNo */
    5,                                                             /* colNo */
    "find",                                                        /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pName */
};

static emlrtRTEInfo ic_emlrtRTEI = {
    253,                                                           /* lineNo */
    5,                                                             /* colNo */
    "find",                                                        /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pName */
};

static emlrtRTEInfo jc_emlrtRTEI = {
    212,                                                           /* lineNo */
    5,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo kc_emlrtRTEI = {
    216,                                                           /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo lc_emlrtRTEI = {
    222,                                                           /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo mc_emlrtRTEI = {
    60,                                                             /* lineNo */
    20,                                                             /* colNo */
    "dot",                                                          /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/specfun/dot.m" /* pName */
};

static emlrtRTEInfo nc_emlrtRTEI = {
    227,                                                           /* lineNo */
    16,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo oc_emlrtRTEI = {
    227,                                                           /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo pc_emlrtRTEI = {
    223,                                                           /* lineNo */
    1,                                                             /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo qc_emlrtRTEI = {
    192,                                                           /* lineNo */
    25,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo rc_emlrtRTEI = {
    247,                                                           /* lineNo */
    17,                                                            /* colNo */
    "find",                                                        /* fName */
    "/usr/local/MATLAB/R2024a/toolbox/eml/lib/matlab/elmat/find.m" /* pName */
};

static emlrtRTEInfo sc_emlrtRTEI = {
    223,                                                           /* lineNo */
    19,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

static emlrtRTEInfo dd_emlrtRTEI = {
    91,                                                            /* lineNo */
    29,                                                            /* colNo */
    "generateAstralNetwork",                                       /* fName */
    "/home/bcberg/Documents/astral-mikado/generateAstralNetwork.m" /* pName */
};

/* Function Declarations */
static void binary_expand_op(const emlrtStack *sp, emxArray_real_T *in1,
                             const emlrtRSInfo in2, const emxArray_real_T *in4,
                             const emxArray_real_T *in5);

static void findNodes(const emlrtStack *sp, const emxArray_real_T *centers,
                      const emxArray_real_T *orients, real_T l,
                      emxArray_real_T *nodes, emxArray_real_T *filCross);

static void sortNodes(const emlrtStack *sp, real_T filIdx,
                      const emxArray_real_T *nodes,
                      const emxArray_real_T *filCross,
                      const emxArray_real_T *centers, real_T astralNum,
                      emxArray_real_T *nodeOrdering);

/* Function Definitions */
static void binary_expand_op(const emlrtStack *sp, emxArray_real_T *in1,
                             const emlrtRSInfo in2, const emxArray_real_T *in4,
                             const emxArray_real_T *in5)
{
  jmp_buf *volatile emlrtJBStack;
  emlrtStack st;
  emxArray_real_T *b_in4;
  emxArray_real_T *r;
  const real_T *in4_data;
  const real_T *in5_data;
  real_T b_varargin_1;
  real_T *b_in4_data;
  real_T *r1;
  int32_T b_loop_ub;
  int32_T binary_expand_op_numThreads;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T loop_ub;
  int32_T stride_0_0;
  int32_T stride_1_0;
  st.prev = sp;
  st.tls = sp->tls;
  in5_data = in5->data;
  in4_data = in4->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in4, 2, &dd_emlrtRTEI);
  if (in5->size[0] == 1) {
    loop_ub = in4->size[0];
  } else {
    loop_ub = in5->size[0];
  }
  i = b_in4->size[0] * b_in4->size[1];
  b_in4->size[0] = loop_ub;
  b_in4->size[1] = 2;
  emxEnsureCapacity_real_T(sp, b_in4, i, &dd_emlrtRTEI);
  b_in4_data = b_in4->data;
  stride_0_0 = (in4->size[0] != 1);
  stride_1_0 = (in5->size[0] != 1);
  b_loop_ub = loop_ub << 1;
  i = (b_loop_ub < 3200);
  if (i) {
    for (i1 = 0; i1 < 2; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_in4_data[i2 + b_in4->size[0] * i1] =
            in4_data[i2 * stride_0_0 + in4->size[0] * i1] -
            in5_data[i2 * stride_1_0 + in5->size[0] * i1];
      }
    }
  } else {
    emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
    emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    binary_expand_op_numThreads = emlrtAllocRegionTLSs(
        sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(binary_expand_op_numThreads) private(i2)

    for (i1 = 0; i1 < 2; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_in4_data[i2 + b_in4->size[0] * i1] =
            in4_data[i2 * stride_0_0 + in4->size[0] * i1] -
            in5_data[i2 * stride_1_0 + in5->size[0] * i1];
      }
    }
    emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  }
  emxInit_real_T(sp, &r, 2, &nb_emlrtRTEI);
  stride_0_0 = r->size[0] * r->size[1];
  r->size[0] = loop_ub;
  r->size[1] = 2;
  emxEnsureCapacity_real_T(sp, r, stride_0_0, &nb_emlrtRTEI);
  r1 = r->data;
  if (i) {
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      real_T varargin_1;
      varargin_1 = b_in4_data[i1];
      r1[i1] = varargin_1 * varargin_1;
    }
  } else {
    emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
    emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    binary_expand_op_numThreads = emlrtAllocRegionTLSs(
        sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(binary_expand_op_numThreads) private(     \
        b_varargin_1)

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_varargin_1 = b_in4_data[i1];
      r1[i1] = b_varargin_1 * b_varargin_1;
    }
    emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  }
  emxFree_real_T(sp, &b_in4);
  st.site = (emlrtRSInfo *)&in2;
  sum(&st, r, in1);
  emxFree_real_T(sp, &r);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void findNodes(const emlrtStack *sp, const emxArray_real_T *centers,
                      const emxArray_real_T *orients, real_T l,
                      emxArray_real_T *nodes, emxArray_real_T *filCross)
{
  jmp_buf *volatile emlrtJBStack;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  emxArray_boolean_T *tooFar;
  emxArray_int32_T *r2;
  emxArray_real_T *centerSepSQR;
  emxArray_real_T *cosines;
  emxArray_real_T *oldNodes;
  emxArray_real_T *otherAsters;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  emxArray_real_T *sines;
  real_T b_otherAsters[2];
  const real_T *centers_data;
  const real_T *orients_data;
  real_T a;
  real_T c;
  real_T numNodesGuess;
  real_T *cosines_data;
  real_T *filCross_data;
  real_T *nodes_data;
  real_T *oldNodes_data;
  real_T *otherAsters_data;
  real_T *r3;
  real_T *r4;
  real_T *sines_data;
  int32_T iv[2];
  int32_T astralNum;
  int32_T astralNum_tmp;
  int32_T b_centers;
  int32_T findNodes_numThreads;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  int32_T idx;
  int32_T jdx;
  int32_T k;
  int32_T loop_ub;
  int32_T loop_ub_tmp;
  int32_T numAsters;
  int32_T numAsters_tmp;
  int32_T numFil;
  int32_T *r5;
  uint32_T nodeCount;
  boolean_T *tooFar_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  orients_data = orients->data;
  centers_data = centers->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /*  Auxiliary functions */
  /*  FINDNODES identifies where filaments cross and records which filaments */
  /*  crossed  */
  /*  (Note: filaments are indexed as one would read left to right across the */
  /*  rows of 'orients', and each row represents an aster. This function does */
  /*  NOT search for intersections between filaments on the same aster) */
  /*    Inputs: */
  /*        centers (numAsters x 2 double): (x,y) coordinates of astral centers
   */
  /*        orients (numAsters x astralNum double): angles giving each  */
  /*        filament's orientation about its astral center */
  /*        l (scalar): length of individual filament */
  /*        centerCross (numAsters x astralNum double): filament indices */
  /*        grouped by each aster */
  /*    Outputs: */
  /*        nodes (numNodes x 2 double): list of (x,y) coordinates of filament
   */
  /*        crossings (EXCLUDING astral centers) */
  /*        filCross (numNodes x 2 double): list of pairs of filament indices */
  /*        corresponding to the filaments that cross at a particular node; */
  /*        listed so that filCross(idx,1) < filCross(idx,2) */
  astralNum_tmp = orients->size[1];
  astralNum = orients->size[1];
  numAsters_tmp = orients->size[0];
  numAsters = orients->size[0];
  numFil = orients->size[0] * orients->size[1];
  st.site = &q_emlrtRSI;
  emxInit_real_T(&st, &cosines, 2, &db_emlrtRTEI);
  i = cosines->size[0] * cosines->size[1];
  cosines->size[0] = orients->size[0];
  cosines->size[1] = orients->size[1];
  emxEnsureCapacity_real_T(&st, cosines, i, &db_emlrtRTEI);
  cosines_data = cosines->data;
  i = (numFil < 3200);
  if (i) {
    for (k = 0; k < numFil; k++) {
      cosines_data[k] = orients_data[k];
    }
  } else {
    emlrtEnterParallelRegion(&st, omp_in_parallel());
    emlrtPushJmpBuf(&st, &emlrtJBStack);
    findNodes_numThreads = emlrtAllocRegionTLSs(
        st.tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(findNodes_numThreads)

    for (k = 0; k < numFil; k++) {
      cosines_data[k] = orients_data[k];
    }
    emlrtPopJmpBuf(&st, &emlrtJBStack);
    emlrtExitParallelRegion(&st, omp_in_parallel());
  }
  b_st.site = &r_emlrtRSI;
  c_st.site = &s_emlrtRSI;
  if (numFil > 2147483646) {
    d_st.site = &t_emlrtRSI;
    check_forloop_overflow_error(&d_st);
  }
  if (i) {
    for (k = 0; k < numFil; k++) {
      cosines_data[k] = muDoubleScalarCos(cosines_data[k]);
    }
  } else {
    emlrtEnterParallelRegion(&b_st, omp_in_parallel());
    emlrtPushJmpBuf(&b_st, &emlrtJBStack);
    findNodes_numThreads =
        emlrtAllocRegionTLSs(b_st.tls, omp_in_parallel(), omp_get_max_threads(),
                             omp_get_num_procs());
#pragma omp parallel for num_threads(findNodes_numThreads)

    for (k = 0; k < numFil; k++) {
      cosines_data[k] = muDoubleScalarCos(cosines_data[k]);
    }
    emlrtPopJmpBuf(&b_st, &emlrtJBStack);
    emlrtExitParallelRegion(&b_st, omp_in_parallel());
  }
  st.site = &p_emlrtRSI;
  emxInit_real_T(&st, &sines, 2, &eb_emlrtRTEI);
  i1 = sines->size[0] * sines->size[1];
  sines->size[0] = numAsters_tmp;
  sines->size[1] = astralNum_tmp;
  emxEnsureCapacity_real_T(&st, sines, i1, &eb_emlrtRTEI);
  sines_data = sines->data;
  if (i) {
    for (k = 0; k < numFil; k++) {
      sines_data[k] = orients_data[k];
    }
  } else {
    emlrtEnterParallelRegion(&st, omp_in_parallel());
    emlrtPushJmpBuf(&st, &emlrtJBStack);
    findNodes_numThreads = emlrtAllocRegionTLSs(
        st.tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(findNodes_numThreads)

    for (k = 0; k < numFil; k++) {
      sines_data[k] = orients_data[k];
    }
    emlrtPopJmpBuf(&st, &emlrtJBStack);
    emlrtExitParallelRegion(&st, omp_in_parallel());
  }
  b_st.site = &u_emlrtRSI;
  c_st.site = &s_emlrtRSI;
  if (i) {
    for (k = 0; k < numFil; k++) {
      sines_data[k] = muDoubleScalarSin(sines_data[k]);
    }
  } else {
    emlrtEnterParallelRegion(&b_st, omp_in_parallel());
    emlrtPushJmpBuf(&b_st, &emlrtJBStack);
    findNodes_numThreads =
        emlrtAllocRegionTLSs(b_st.tls, omp_in_parallel(), omp_get_max_threads(),
                             omp_get_num_procs());
#pragma omp parallel for num_threads(findNodes_numThreads)

    for (k = 0; k < numFil; k++) {
      sines_data[k] = muDoubleScalarSin(sines_data[k]);
    }
    emlrtPopJmpBuf(&b_st, &emlrtJBStack);
    emlrtExitParallelRegion(&b_st, omp_in_parallel());
  }
  numNodesGuess = 4.0 * l * (real_T)numFil;
  /*  Note: this is an empirical guess! */
  if (!(numNodesGuess >= 0.0)) {
    emlrtNonNegativeCheckR2012b(numNodesGuess, &h_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = (int32_T)muDoubleScalarFloor(numNodesGuess);
  if (numNodesGuess != i) {
    emlrtIntegerCheckR2012b(numNodesGuess, &g_emlrtDCI, (emlrtConstCTX)sp);
  }
  i1 = (int32_T)numNodesGuess;
  i2 = nodes->size[0] * nodes->size[1];
  nodes->size[0] = (int32_T)numNodesGuess;
  nodes->size[1] = 2;
  emxEnsureCapacity_real_T(sp, nodes, i2, &fb_emlrtRTEI);
  nodes_data = nodes->data;
  if (numNodesGuess != i) {
    emlrtIntegerCheckR2012b(numNodesGuess, &n_emlrtDCI, (emlrtConstCTX)sp);
  }
  loop_ub_tmp = (int32_T)numNodesGuess << 1;
  if (loop_ub_tmp < 3200) {
    for (k = 0; k < loop_ub_tmp; k++) {
      nodes_data[k] = 0.0;
    }
  } else {
    emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
    emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    findNodes_numThreads = emlrtAllocRegionTLSs(
        sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(findNodes_numThreads)

    for (k = 0; k < loop_ub_tmp; k++) {
      nodes_data[k] = 0.0;
    }
    emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  }
  if ((int32_T)numNodesGuess != i) {
    emlrtIntegerCheckR2012b(numNodesGuess, &i_emlrtDCI, (emlrtConstCTX)sp);
  }
  i2 = filCross->size[0] * filCross->size[1];
  filCross->size[0] = i1;
  filCross->size[1] = 2;
  emxEnsureCapacity_real_T(sp, filCross, i2, &gb_emlrtRTEI);
  filCross_data = filCross->data;
  if ((int32_T)numNodesGuess != i) {
    emlrtIntegerCheckR2012b(numNodesGuess, &o_emlrtDCI, (emlrtConstCTX)sp);
  }
  if (loop_ub_tmp < 3200) {
    for (k = 0; k < loop_ub_tmp; k++) {
      filCross_data[k] = 0.0;
    }
  } else {
    emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
    emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    findNodes_numThreads = emlrtAllocRegionTLSs(
        sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(findNodes_numThreads)

    for (k = 0; k < loop_ub_tmp; k++) {
      filCross_data[k] = 0.0;
    }
    emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  }
  nodeCount = 0U;
  emxInit_real_T(sp, &centerSepSQR, 1, &yb_emlrtRTEI);
  emxInit_real_T(sp, &oldNodes, 2, &qb_emlrtRTEI);
  emxInit_real_T(sp, &otherAsters, 2, &jb_emlrtRTEI);
  otherAsters_data = otherAsters->data;
  emxInit_real_T(sp, &r, 2, &ac_emlrtRTEI);
  emxInit_real_T(sp, &r1, 2, &nb_emlrtRTEI);
  if (orients->size[1] == 1) {
    /*  routine for "Classical Mikado" networks */
    if (numFil - 2 >= 0) {
      b_otherAsters[1] = 1.0;
      b_centers = centers->size[0];
      a = 2.0 * l;
    }
    emxInit_int32_T(sp, &r2, &bc_emlrtRTEI);
    for (idx = 0; idx <= numFil - 2; idx++) {
      real_T c_centers[2];
      int32_T b_loop_ub;
      /*  only look for intersections if filaments are close enough */
      st.site = &o_emlrtRSI;
      b_st.site = &v_emlrtRSI;
      if (numFil < ((real_T)idx + 1.0) + 1.0) {
        otherAsters->size[0] = 1;
        otherAsters->size[1] = 0;
      } else {
        i = otherAsters->size[0] * otherAsters->size[1];
        otherAsters->size[0] = 1;
        i1 = numFil - idx;
        otherAsters->size[1] = i1 - 1;
        emxEnsureCapacity_real_T(&b_st, otherAsters, i, &ib_emlrtRTEI);
        otherAsters_data = otherAsters->data;
        loop_ub = i1 - 2;
        for (i = 0; i <= loop_ub; i++) {
          otherAsters_data[i] = (real_T)((uint32_T)idx + (uint32_T)i) + 2.0;
        }
      }
      i = centers->size[0];
      if (((int32_T)((uint32_T)idx + 1U) < 1) ||
          ((int32_T)((uint32_T)idx + 1U) > centers->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)idx + 1U), 1,
                                      centers->size[0], &nc_emlrtBCI,
                                      (emlrtConstCTX)sp);
      }
      c_centers[0] = centers_data[idx];
      c_centers[1] = centers_data[idx + centers->size[0]];
      b_otherAsters[0] = otherAsters->size[1];
      st.site = &n_emlrtRSI;
      repmat(&st, c_centers, b_otherAsters, oldNodes);
      oldNodes_data = oldNodes->data;
      loop_ub = otherAsters->size[1];
      i1 = r->size[0] * r->size[1];
      r->size[0] = otherAsters->size[1];
      r->size[1] = 2;
      emxEnsureCapacity_real_T(sp, r, i1, &lb_emlrtRTEI);
      r3 = r->data;
      for (i1 = 0; i1 < 2; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          i3 = (int32_T)otherAsters_data[i2];
          if (otherAsters_data[i2] != i3) {
            emlrtIntegerCheckR2012b(otherAsters_data[i2], &p_emlrtDCI,
                                    (emlrtConstCTX)sp);
          }
          if ((i3 < 1) || (i3 > b_centers)) {
            emlrtDynamicBoundsCheckR2012b(i3, 1, b_centers, &oc_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          r3[i2 + r->size[0] * i1] =
              centers_data[(i3 + centers->size[0] * i1) - 1];
        }
      }
      if ((oldNodes->size[0] != otherAsters->size[1]) &&
          ((oldNodes->size[0] != 1) && (otherAsters->size[1] != 1))) {
        emlrtDimSizeImpxCheckR2021b(oldNodes->size[0], otherAsters->size[1],
                                    &g_emlrtECI, (emlrtConstCTX)sp);
      }
      if (oldNodes->size[0] == r->size[0]) {
        i1 = r1->size[0] * r1->size[1];
        r1->size[0] = oldNodes->size[0];
        r1->size[1] = 2;
        emxEnsureCapacity_real_T(sp, r1, i1, &nb_emlrtRTEI);
        r4 = r1->data;
        loop_ub = oldNodes->size[0] << 1;
        for (i1 = 0; i1 < loop_ub; i1++) {
          real_T filSubJdx;
          filSubJdx = oldNodes_data[i1] - r3[i1];
          r4[i1] = filSubJdx * filSubJdx;
        }
        st.site = &n_emlrtRSI;
        sum(&st, r1, centerSepSQR);
        oldNodes_data = centerSepSQR->data;
      } else {
        st.site = &n_emlrtRSI;
        binary_expand_op(&st, centerSepSQR, n_emlrtRSI, oldNodes, r);
        oldNodes_data = centerSepSQR->data;
      }
      st.site = &m_emlrtRSI;
      b_st.site = &hb_emlrtRSI;
      c_st.site = &ib_emlrtRSI;
      loop_ub = centerSepSQR->size[0] - 1;
      b_loop_ub = 0;
      for (astralNum_tmp = 0; astralNum_tmp <= loop_ub; astralNum_tmp++) {
        if (oldNodes_data[astralNum_tmp] <= a * a) {
          b_loop_ub++;
        }
      }
      i1 = r2->size[0];
      r2->size[0] = b_loop_ub;
      emxEnsureCapacity_int32_T(sp, r2, i1, &ob_emlrtRTEI);
      r5 = r2->data;
      b_loop_ub = 0;
      for (astralNum_tmp = 0; astralNum_tmp <= loop_ub; astralNum_tmp++) {
        if (oldNodes_data[astralNum_tmp] <= a * a) {
          r5[b_loop_ub] = astralNum_tmp;
          b_loop_ub++;
        }
      }
      loop_ub_tmp = r2->size[0];
      for (i1 = 0; i1 < loop_ub_tmp; i1++) {
        if (r5[i1] > otherAsters->size[1] - 1) {
          emlrtDynamicBoundsCheckR2012b(r5[i1], 0, otherAsters->size[1] - 1,
                                        &qc_emlrtBCI, (emlrtConstCTX)sp);
        }
      }
      for (jdx = 0; jdx < loop_ub_tmp; jdx++) {
        real_T b1;
        real_T b2;
        real_T denom;
        real_T denom_tmp;
        real_T t1;
        if (jdx + 1 > loop_ub_tmp) {
          emlrtDynamicBoundsCheckR2012b(jdx + 1, 1, loop_ub_tmp, &mc_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        /*  legacy code in these comments for readability */
        /*  A = [cosines(idx), - cosines(filOfJdx); */
        /*      sines(idx), - sines(filOfJdx)]; */
        /*  denom = det(A); */
        /*  b = 1/l * [centers(filOfJdx,1) - centers(idx,1); */
        /*      centers(filOfJdx,2) - centers(idx,2)]; */
        /*  t1 = det([b,A(:,2)]) / denom; */
        /*  t2 = det([A(:,1),b]) / denom; */
        if (((int32_T)((uint32_T)idx + 1U) < 1) ||
            ((int32_T)((uint32_T)idx + 1U) > numFil)) {
          emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)idx + 1U), 1,
                                        numFil, &sc_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        i1 = (int32_T)otherAsters_data[r5[jdx]];
        if ((i1 < 1) || (i1 > numFil)) {
          emlrtDynamicBoundsCheckR2012b(i1, 1, numFil, &tc_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        if (((int32_T)((uint32_T)idx + 1U) < 1) ||
            ((int32_T)((uint32_T)idx + 1U) > numFil)) {
          emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)idx + 1U), 1,
                                        numFil, &uc_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        if (i1 > numFil) {
          emlrtDynamicBoundsCheckR2012b(i1, 1, numFil, &vc_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        t1 = -sines_data[i1 - 1];
        denom_tmp = cosines_data[i1 - 1];
        denom = cosines_data[idx] * t1 + sines_data[idx] * denom_tmp;
        if (i1 > i) {
          emlrtDynamicBoundsCheckR2012b(i1, 1, i, &wc_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        if (((int32_T)((uint32_T)idx + 1U) < 1) ||
            ((int32_T)((uint32_T)idx + 1U) > i)) {
          emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)idx + 1U), 1, i,
                                        &xc_emlrtBCI, (emlrtConstCTX)sp);
        }
        b1 = (centers_data[i1 - 1] - centers_data[idx]) / l;
        if (((int32_T)((uint32_T)idx + 1U) < 1) ||
            ((int32_T)((uint32_T)idx + 1U) > i)) {
          emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)idx + 1U), 1, i,
                                        &yc_emlrtBCI, (emlrtConstCTX)sp);
        }
        b2 = (centers_data[(i1 + centers->size[0]) - 1] -
              centers_data[idx + centers->size[0]]) /
             l;
        if (i1 > numFil) {
          emlrtDynamicBoundsCheckR2012b(i1, 1, numFil, &ad_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        t1 = (b1 * t1 + b2 * denom_tmp) / denom;
        if ((idx + 1 < 1) || (idx + 1 > numFil)) {
          emlrtDynamicBoundsCheckR2012b(idx + 1, 1, numFil, &lc_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        if ((muDoubleScalarAbs(t1 - 0.5) <= 0.5) &&
            (muDoubleScalarAbs((cosines_data[idx] * b2 - sines_data[idx] * b1) /
                                   denom -
                               0.5) <= 0.5)) {
          nodeCount++;
          if (nodeCount > numNodesGuess) {
            numNodesGuess += (real_T)numFil;
            i1 = oldNodes->size[0] * oldNodes->size[1];
            oldNodes->size[0] = nodes->size[0];
            oldNodes->size[1] = 2;
            emxEnsureCapacity_real_T(sp, oldNodes, i1, &qb_emlrtRTEI);
            oldNodes_data = oldNodes->data;
            loop_ub = nodes->size[0] << 1;
            for (i1 = 0; i1 < loop_ub; i1++) {
              oldNodes_data[i1] = nodes_data[i1];
            }
            if (numNodesGuess != (int32_T)numNodesGuess) {
              emlrtIntegerCheckR2012b(numNodesGuess, &j_emlrtDCI,
                                      (emlrtConstCTX)sp);
            }
            i1 = nodes->size[0] * nodes->size[1];
            nodes->size[0] = (int32_T)numNodesGuess;
            nodes->size[1] = 2;
            emxEnsureCapacity_real_T(sp, nodes, i1, &rb_emlrtRTEI);
            nodes_data = nodes->data;
            if ((int32_T)numNodesGuess != (int32_T)numNodesGuess) {
              emlrtIntegerCheckR2012b(numNodesGuess, &t_emlrtDCI,
                                      (emlrtConstCTX)sp);
            }
            b_loop_ub = (int32_T)numNodesGuess << 1;
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              nodes_data[i1] = 0.0;
            }
            if ((int32_T)(nodeCount - 1U) < 1) {
              loop_ub = 0;
            } else {
              if (numNodesGuess < 1.0) {
                emlrtDynamicBoundsCheckR2012b(1, 1, (int32_T)numNodesGuess,
                                              &kc_emlrtBCI, (emlrtConstCTX)sp);
              }
              if (((int32_T)(nodeCount - 1U) < 1) ||
                  ((int32_T)(nodeCount - 1U) > (int32_T)numNodesGuess)) {
                emlrtDynamicBoundsCheckR2012b((int32_T)(nodeCount - 1U), 1,
                                              (int32_T)numNodesGuess,
                                              &jc_emlrtBCI, (emlrtConstCTX)sp);
              }
              loop_ub = (int32_T)(nodeCount - 1U);
            }
            iv[0] = loop_ub;
            iv[1] = 2;
            emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &oldNodes->size[0], 2,
                                          &f_emlrtECI, (emlrtCTX)sp);
            for (i1 = 0; i1 < 2; i1++) {
              for (i2 = 0; i2 < loop_ub; i2++) {
                nodes_data[i2 + nodes->size[0] * i1] =
                    oldNodes_data[i2 + oldNodes->size[0] * i1];
              }
            }
            i1 = oldNodes->size[0] * oldNodes->size[1];
            oldNodes->size[0] = filCross->size[0];
            oldNodes->size[1] = 2;
            emxEnsureCapacity_real_T(sp, oldNodes, i1, &sb_emlrtRTEI);
            oldNodes_data = oldNodes->data;
            loop_ub = filCross->size[0] << 1;
            for (i1 = 0; i1 < loop_ub; i1++) {
              oldNodes_data[i1] = filCross_data[i1];
            }
            i1 = (int32_T)numNodesGuess;
            if (i1 != (int32_T)numNodesGuess) {
              emlrtIntegerCheckR2012b(numNodesGuess, &k_emlrtDCI,
                                      (emlrtConstCTX)sp);
            }
            i2 = filCross->size[0] * filCross->size[1];
            filCross->size[0] = i1;
            filCross->size[1] = 2;
            emxEnsureCapacity_real_T(sp, filCross, i2, &tb_emlrtRTEI);
            filCross_data = filCross->data;
            if (i1 != (int32_T)numNodesGuess) {
              emlrtIntegerCheckR2012b(numNodesGuess, &v_emlrtDCI,
                                      (emlrtConstCTX)sp);
            }
            for (i2 = 0; i2 < b_loop_ub; i2++) {
              filCross_data[i2] = 0.0;
            }
            if ((int32_T)(nodeCount - 1U) < 1) {
              loop_ub = 0;
            } else {
              if (numNodesGuess < 1.0) {
                emlrtDynamicBoundsCheckR2012b(1, 1, (int32_T)numNodesGuess,
                                              &ic_emlrtBCI, (emlrtConstCTX)sp);
              }
              if (((int32_T)(nodeCount - 1U) < 1) ||
                  ((int32_T)(nodeCount - 1U) > i1)) {
                emlrtDynamicBoundsCheckR2012b((int32_T)(nodeCount - 1U), 1,
                                              (int32_T)numNodesGuess,
                                              &hc_emlrtBCI, (emlrtConstCTX)sp);
              }
              loop_ub = (int32_T)(nodeCount - 1U);
            }
            iv[0] = loop_ub;
            iv[1] = 2;
            emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &oldNodes->size[0], 2,
                                          &e_emlrtECI, (emlrtCTX)sp);
            for (i1 = 0; i1 < 2; i1++) {
              for (i2 = 0; i2 < loop_ub; i2++) {
                filCross_data[i2 + filCross->size[0] * i1] =
                    oldNodes_data[i2 + oldNodes->size[0] * i1];
              }
            }
          }
          if ((idx + 1 < 1) || (idx + 1 > i)) {
            emlrtDynamicBoundsCheckR2012b(idx + 1, 1, i, &gc_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          if ((idx + 1 < 1) || (idx + 1 > numFil)) {
            emlrtDynamicBoundsCheckR2012b(idx + 1, 1, numFil, &fc_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          if ((idx + 1 < 1) || (idx + 1 > i)) {
            emlrtDynamicBoundsCheckR2012b(idx + 1, 1, i, &ec_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          if ((idx + 1 < 1) || (idx + 1 > numFil)) {
            emlrtDynamicBoundsCheckR2012b(idx + 1, 1, numFil, &dc_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          if (((int32_T)nodeCount < 1) ||
              ((int32_T)nodeCount > nodes->size[0])) {
            emlrtDynamicBoundsCheckR2012b((int32_T)nodeCount, 1, nodes->size[0],
                                          &cc_emlrtBCI, (emlrtConstCTX)sp);
          }
          nodes_data[(int32_T)nodeCount - 1] =
              centers_data[idx] + l * cosines_data[idx] * t1;
          nodes_data[((int32_T)nodeCount + nodes->size[0]) - 1] =
              centers_data[idx + centers->size[0]] + l * sines_data[idx] * t1;
          if (((int32_T)nodeCount < 1) ||
              ((int32_T)nodeCount > filCross->size[0])) {
            emlrtDynamicBoundsCheckR2012b((int32_T)nodeCount, 1,
                                          filCross->size[0], &bc_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          filCross_data[(int32_T)nodeCount - 1] = (real_T)idx + 1.0;
          filCross_data[((int32_T)nodeCount + filCross->size[0]) - 1] =
              otherAsters_data[r5[jdx]];
        }
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b((emlrtConstCTX)sp);
        }
      }
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b((emlrtConstCTX)sp);
      }
    }
    emxFree_int32_T(sp, &r2);
  } else if (orients->size[1] >= 2) {
    real_T d;
    /*  routine for "Astral Mikado" networks */
    /*  loop only needs to run through the last filament on penultimate aster */
    d = ((real_T)orients->size[0] - 1.0) * (real_T)orients->size[1];
    i = (int32_T)d;
    emlrtForLoopVectorCheckR2021a(1.0, 1.0, d, mxDOUBLE_CLASS, (int32_T)d,
                                  &h_emlrtRTEI, (emlrtConstCTX)sp);
    if ((int32_T)d - 1 >= 0) {
      b_otherAsters[1] = 1.0;
      b_centers = centers->size[0];
      a = 2.0 * l;
      c = a * a;
    }
    emxInit_boolean_T(sp, &tooFar, 1, &pb_emlrtRTEI);
    for (idx = 0; idx < i; idx++) {
      real_T c_centers[2];
      real_T b_centers_tmp;
      real_T centers_tmp;
      real_T filSubJdx;
      real_T y;
      boolean_T b;
      y = muDoubleScalarFloor((((real_T)idx + 1.0) - 1.0) / (real_T)astralNum);
      st.site = &l_emlrtRSI;
      b_st.site = &v_emlrtRSI;
      if (numAsters < (y + 1.0) + 1.0) {
        otherAsters->size[0] = 1;
        otherAsters->size[1] = 0;
      } else if ((y + 1.0) + 1.0 == (y + 1.0) + 1.0) {
        i1 = otherAsters->size[0] * otherAsters->size[1];
        otherAsters->size[0] = 1;
        loop_ub = (int32_T)((real_T)numAsters - ((y + 1.0) + 1.0));
        otherAsters->size[1] = loop_ub + 1;
        emxEnsureCapacity_real_T(&b_st, otherAsters, i1, &jb_emlrtRTEI);
        otherAsters_data = otherAsters->data;
        for (i1 = 0; i1 <= loop_ub; i1++) {
          otherAsters_data[i1] = ((y + 1.0) + 1.0) + (real_T)i1;
        }
      } else {
        c_st.site = &w_emlrtRSI;
        eml_float_colon(&c_st, (y + 1.0) + 1.0, numAsters, otherAsters);
        otherAsters_data = otherAsters->data;
      }
      i1 = centers->size[0];
      b = (((int32_T)(y + 1.0) < 1) || ((int32_T)(y + 1.0) > centers->size[0]));
      if (b) {
        emlrtDynamicBoundsCheckR2012b((int32_T)(y + 1.0), 1, centers->size[0],
                                      &ac_emlrtBCI, (emlrtConstCTX)sp);
      }
      centers_tmp = centers_data[(int32_T)(y + 1.0) - 1];
      c_centers[0] = centers_tmp;
      b_centers_tmp = centers_data[((int32_T)(y + 1.0) + centers->size[0]) - 1];
      c_centers[1] = b_centers_tmp;
      b_otherAsters[0] = otherAsters->size[1];
      st.site = &k_emlrtRSI;
      repmat(&st, c_centers, b_otherAsters, oldNodes);
      oldNodes_data = oldNodes->data;
      loop_ub = otherAsters->size[1];
      i2 = r->size[0] * r->size[1];
      r->size[0] = otherAsters->size[1];
      r->size[1] = 2;
      emxEnsureCapacity_real_T(sp, r, i2, &mb_emlrtRTEI);
      r3 = r->data;
      for (i2 = 0; i2 < 2; i2++) {
        for (i3 = 0; i3 < loop_ub; i3++) {
          if (otherAsters_data[i3] !=
              (int32_T)muDoubleScalarFloor(otherAsters_data[i3])) {
            emlrtIntegerCheckR2012b(otherAsters_data[i3], &q_emlrtDCI,
                                    (emlrtConstCTX)sp);
          }
          i4 = (int32_T)otherAsters_data[i3];
          if ((i4 < 1) || (i4 > b_centers)) {
            emlrtDynamicBoundsCheckR2012b(i4, 1, b_centers, &pc_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          r3[i3 + r->size[0] * i2] =
              centers_data[(i4 + centers->size[0] * i2) - 1];
        }
      }
      if ((oldNodes->size[0] != otherAsters->size[1]) &&
          ((oldNodes->size[0] != 1) && (otherAsters->size[1] != 1))) {
        emlrtDimSizeImpxCheckR2021b(oldNodes->size[0], otherAsters->size[1],
                                    &d_emlrtECI, (emlrtConstCTX)sp);
      }
      if (oldNodes->size[0] == r->size[0]) {
        i2 = r1->size[0] * r1->size[1];
        r1->size[0] = oldNodes->size[0];
        r1->size[1] = 2;
        emxEnsureCapacity_real_T(sp, r1, i2, &nb_emlrtRTEI);
        r4 = r1->data;
        loop_ub = oldNodes->size[0] << 1;
        for (i2 = 0; i2 < loop_ub; i2++) {
          filSubJdx = oldNodes_data[i2] - r3[i2];
          r4[i2] = filSubJdx * filSubJdx;
        }
        st.site = &k_emlrtRSI;
        sum(&st, r1, centerSepSQR);
        oldNodes_data = centerSepSQR->data;
      } else {
        st.site = &k_emlrtRSI;
        binary_expand_op(&st, centerSepSQR, k_emlrtRSI, oldNodes, r);
        oldNodes_data = centerSepSQR->data;
      }
      st.site = &j_emlrtRSI;
      b_st.site = &hb_emlrtRSI;
      c_st.site = &ib_emlrtRSI;
      loop_ub = centerSepSQR->size[0];
      i2 = tooFar->size[0];
      tooFar->size[0] = centerSepSQR->size[0];
      emxEnsureCapacity_boolean_T(sp, tooFar, i2, &pb_emlrtRTEI);
      tooFar_data = tooFar->data;
      for (i2 = 0; i2 < loop_ub; i2++) {
        tooFar_data[i2] = (oldNodes_data[i2] > c);
      }
      /*  start second loop at first filament on the next aster */
      d = (y + 1.0) * (real_T)astralNum + 1.0;
      i2 = (int32_T)((real_T)numFil + (1.0 - d));
      emlrtForLoopVectorCheckR2021a(d, 1.0, numFil, mxDOUBLE_CLASS, i2,
                                    &g_emlrtRTEI, (emlrtConstCTX)sp);
      for (jdx = 0; jdx < i2; jdx++) {
        real_T b_jdx;
        b_jdx = d + (real_T)jdx;
        a = muDoubleScalarFloor((b_jdx - 1.0) / (real_T)astralNum);
        /*  whichOtherAster finds correct row in centerSepSQR, tooFar */
        i3 = (int32_T)((a + 1.0) - (y + 1.0));
        if ((i3 < 1) || (i3 > loop_ub)) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, loop_ub, &rc_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        if (!tooFar_data[i3 - 1]) {
          real_T b1;
          real_T b2;
          real_T denom;
          real_T denom_tmp;
          real_T filSubIdx;
          real_T t1;
          /*  otherwise, check for intersection */
          filSubJdx = b_mod((real_T)idx + 1.0, astralNum);
          filSubIdx = filSubJdx + (real_T)(astralNum * (filSubJdx == 0.0));
          filSubJdx = b_mod(b_jdx, astralNum);
          filSubJdx += (real_T)(astralNum * (filSubJdx == 0.0));
          /*  legacy code in these comments for readability */
          /*  A = [cosines(asterIdx,filSubIdx), - cosines(asterJdx,filSubJdx);
           */
          /*      sines(asterIdx,filSubIdx), - sines(asterJdx,filSubJdx)]; */
          /*  denom = det(A); */
          /*  b = 1/l * [centers(asterJdx,1) - centers(asterIdx,1); */
          /*      centers(asterJdx,2) - centers(asterIdx,2)]; */
          /*  t1 = det([b,A(:,2)]) / denom; */
          /*  t2 = det([A(:,1),b]) / denom; */
          if (((int32_T)(y + 1.0) < 1) ||
              ((int32_T)(y + 1.0) > numAsters_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(y + 1.0), 1, numAsters_tmp,
                                          &bd_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (filSubIdx != (int32_T)muDoubleScalarFloor(filSubIdx)) {
            emlrtIntegerCheckR2012b(filSubIdx, &r_emlrtDCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)filSubIdx < 1) ||
              ((int32_T)filSubIdx > astralNum_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)filSubIdx, 1, astralNum_tmp,
                                          &cd_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)(a + 1.0) < 1) ||
              ((int32_T)(a + 1.0) > numAsters_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(a + 1.0), 1, numAsters_tmp,
                                          &dd_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (filSubJdx != (int32_T)muDoubleScalarFloor(filSubJdx)) {
            emlrtIntegerCheckR2012b(filSubJdx, &s_emlrtDCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)filSubJdx < 1) ||
              ((int32_T)filSubJdx > astralNum_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)filSubJdx, 1, astralNum_tmp,
                                          &ed_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)(y + 1.0) < 1) ||
              ((int32_T)(y + 1.0) > numAsters_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(y + 1.0), 1, numAsters_tmp,
                                          &fd_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)filSubIdx < 1) ||
              ((int32_T)filSubIdx > astralNum_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)filSubIdx, 1, astralNum_tmp,
                                          &gd_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)(a + 1.0) < 1) ||
              ((int32_T)(a + 1.0) > numAsters_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(a + 1.0), 1, numAsters_tmp,
                                          &hd_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)filSubJdx < 1) ||
              ((int32_T)filSubJdx > astralNum_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)filSubJdx, 1, astralNum_tmp,
                                          &id_emlrtBCI, (emlrtConstCTX)sp);
          }
          t1 = cosines_data[((int32_T)(a + 1.0) +
                             cosines->size[0] * ((int32_T)filSubJdx - 1)) -
                            1];
          denom_tmp = sines_data[((int32_T)(y + 1.0) +
                                  sines->size[0] * ((int32_T)filSubIdx - 1)) -
                                 1];
          denom = cosines_data[((int32_T)(y + 1.0) +
                                cosines->size[0] * ((int32_T)filSubIdx - 1)) -
                               1] *
                      -sines_data[((int32_T)(a + 1.0) +
                                   sines->size[0] * ((int32_T)filSubJdx - 1)) -
                                  1] +
                  denom_tmp * t1;
          if (((int32_T)(a + 1.0) < 1) || ((int32_T)(a + 1.0) > i1)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(a + 1.0), 1, i1,
                                          &jd_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)(y + 1.0) < 1) || ((int32_T)(y + 1.0) > i1)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(y + 1.0), 1, i1,
                                          &kd_emlrtBCI, (emlrtConstCTX)sp);
          }
          b1 = (centers_data[(int32_T)(a + 1.0) - 1] - centers_tmp) / l;
          if (((int32_T)(y + 1.0) < 1) || ((int32_T)(y + 1.0) > i1)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(y + 1.0), 1, i1,
                                          &ld_emlrtBCI, (emlrtConstCTX)sp);
          }
          b2 = (centers_data[((int32_T)(a + 1.0) + centers->size[0]) - 1] -
                b_centers_tmp) /
               l;
          if (((int32_T)(a + 1.0) < 1) ||
              ((int32_T)(a + 1.0) > numAsters_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(a + 1.0), 1, numAsters_tmp,
                                          &md_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)filSubJdx < 1) ||
              ((int32_T)filSubJdx > astralNum_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)filSubJdx, 1, astralNum_tmp,
                                          &nd_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)(a + 1.0) < 1) ||
              ((int32_T)(a + 1.0) > numAsters_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(a + 1.0), 1, numAsters_tmp,
                                          &od_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)filSubJdx < 1) ||
              ((int32_T)filSubJdx > astralNum_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)filSubJdx, 1, astralNum_tmp,
                                          &pd_emlrtBCI, (emlrtConstCTX)sp);
          }
          t1 = (b1 * -sines_data[((int32_T)(a + 1.0) +
                                  sines->size[0] * ((int32_T)filSubJdx - 1)) -
                                 1] +
                b2 * t1) /
               denom;
          if (((int32_T)filSubIdx < 1) ||
              ((int32_T)filSubIdx > astralNum_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)filSubIdx, 1, astralNum_tmp,
                                          &yb_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)(y + 1.0) < 1) ||
              ((int32_T)(y + 1.0) > numAsters_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)(y + 1.0), 1, numAsters_tmp,
                                          &xb_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)filSubIdx < 1) ||
              ((int32_T)filSubIdx > astralNum_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)filSubIdx, 1, astralNum_tmp,
                                          &wb_emlrtBCI, (emlrtConstCTX)sp);
          }
          if (muDoubleScalarAbs(t1 - 0.5) <= 0.5) {
            filSubJdx =
                cosines_data[((int32_T)(y + 1.0) +
                              cosines->size[0] * ((int32_T)filSubIdx - 1)) -
                             1];
            if (muDoubleScalarAbs((filSubJdx * b2 - denom_tmp * b1) / denom -
                                  0.5) <= 0.5) {
              nodeCount++;
              if (nodeCount > numNodesGuess) {
                int32_T b_loop_ub;
                numNodesGuess += (real_T)numFil;
                i3 = oldNodes->size[0] * oldNodes->size[1];
                oldNodes->size[0] = nodes->size[0];
                oldNodes->size[1] = 2;
                emxEnsureCapacity_real_T(sp, oldNodes, i3, &ub_emlrtRTEI);
                oldNodes_data = oldNodes->data;
                b_loop_ub = nodes->size[0] << 1;
                for (i3 = 0; i3 < b_loop_ub; i3++) {
                  oldNodes_data[i3] = nodes_data[i3];
                }
                if (numNodesGuess != (int32_T)numNodesGuess) {
                  emlrtIntegerCheckR2012b(numNodesGuess, &l_emlrtDCI,
                                          (emlrtConstCTX)sp);
                }
                i3 = nodes->size[0] * nodes->size[1];
                nodes->size[0] = (int32_T)numNodesGuess;
                nodes->size[1] = 2;
                emxEnsureCapacity_real_T(sp, nodes, i3, &vb_emlrtRTEI);
                nodes_data = nodes->data;
                if ((int32_T)numNodesGuess != (int32_T)numNodesGuess) {
                  emlrtIntegerCheckR2012b(numNodesGuess, &u_emlrtDCI,
                                          (emlrtConstCTX)sp);
                }
                loop_ub_tmp = (int32_T)numNodesGuess << 1;
                for (i3 = 0; i3 < loop_ub_tmp; i3++) {
                  nodes_data[i3] = 0.0;
                }
                if ((int32_T)(nodeCount - 1U) < 1) {
                  b_loop_ub = 0;
                } else {
                  if (numNodesGuess < 1.0) {
                    emlrtDynamicBoundsCheckR2012b(1, 1, (int32_T)numNodesGuess,
                                                  &vb_emlrtBCI,
                                                  (emlrtConstCTX)sp);
                  }
                  if (((int32_T)(nodeCount - 1U) < 1) ||
                      ((int32_T)(nodeCount - 1U) > (int32_T)numNodesGuess)) {
                    emlrtDynamicBoundsCheckR2012b(
                        (int32_T)(nodeCount - 1U), 1, (int32_T)numNodesGuess,
                        &ub_emlrtBCI, (emlrtConstCTX)sp);
                  }
                  b_loop_ub = (int32_T)(nodeCount - 1U);
                }
                iv[0] = b_loop_ub;
                iv[1] = 2;
                emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &oldNodes->size[0], 2,
                                              &c_emlrtECI, (emlrtCTX)sp);
                for (i3 = 0; i3 < 2; i3++) {
                  for (i4 = 0; i4 < b_loop_ub; i4++) {
                    nodes_data[i4 + nodes->size[0] * i3] =
                        oldNodes_data[i4 + oldNodes->size[0] * i3];
                  }
                }
                i3 = oldNodes->size[0] * oldNodes->size[1];
                oldNodes->size[0] = filCross->size[0];
                oldNodes->size[1] = 2;
                emxEnsureCapacity_real_T(sp, oldNodes, i3, &wb_emlrtRTEI);
                oldNodes_data = oldNodes->data;
                b_loop_ub = filCross->size[0] << 1;
                for (i3 = 0; i3 < b_loop_ub; i3++) {
                  oldNodes_data[i3] = filCross_data[i3];
                }
                if ((int32_T)numNodesGuess != (int32_T)numNodesGuess) {
                  emlrtIntegerCheckR2012b(numNodesGuess, &m_emlrtDCI,
                                          (emlrtConstCTX)sp);
                }
                i3 = filCross->size[0] * filCross->size[1];
                filCross->size[0] = (int32_T)numNodesGuess;
                filCross->size[1] = 2;
                emxEnsureCapacity_real_T(sp, filCross, i3, &xb_emlrtRTEI);
                filCross_data = filCross->data;
                if ((int32_T)numNodesGuess != (int32_T)numNodesGuess) {
                  emlrtIntegerCheckR2012b(numNodesGuess, &w_emlrtDCI,
                                          (emlrtConstCTX)sp);
                }
                for (i3 = 0; i3 < loop_ub_tmp; i3++) {
                  filCross_data[i3] = 0.0;
                }
                if ((int32_T)(nodeCount - 1U) < 1) {
                  b_loop_ub = 0;
                } else {
                  if (numNodesGuess < 1.0) {
                    emlrtDynamicBoundsCheckR2012b(1, 1, (int32_T)numNodesGuess,
                                                  &tb_emlrtBCI,
                                                  (emlrtConstCTX)sp);
                  }
                  if (((int32_T)(nodeCount - 1U) < 1) ||
                      ((int32_T)(nodeCount - 1U) > (int32_T)numNodesGuess)) {
                    emlrtDynamicBoundsCheckR2012b(
                        (int32_T)(nodeCount - 1U), 1, (int32_T)numNodesGuess,
                        &sb_emlrtBCI, (emlrtConstCTX)sp);
                  }
                  b_loop_ub = (int32_T)(nodeCount - 1U);
                }
                iv[0] = b_loop_ub;
                iv[1] = 2;
                emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &oldNodes->size[0], 2,
                                              &b_emlrtECI, (emlrtCTX)sp);
                for (i3 = 0; i3 < 2; i3++) {
                  for (i4 = 0; i4 < b_loop_ub; i4++) {
                    filCross_data[i4 + filCross->size[0] * i3] =
                        oldNodes_data[i4 + oldNodes->size[0] * i3];
                  }
                }
              }
              if (((int32_T)(y + 1.0) < 1) ||
                  ((int32_T)(y + 1.0) > numAsters_tmp)) {
                emlrtDynamicBoundsCheckR2012b((int32_T)(y + 1.0), 1,
                                              numAsters_tmp, &rb_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              if (((int32_T)filSubIdx < 1) ||
                  ((int32_T)filSubIdx > astralNum_tmp)) {
                emlrtDynamicBoundsCheckR2012b((int32_T)filSubIdx, 1,
                                              astralNum_tmp, &qb_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              if (((int32_T)(y + 1.0) < 1) || ((int32_T)(y + 1.0) > i1)) {
                emlrtDynamicBoundsCheckR2012b((int32_T)(y + 1.0), 1, i1,
                                              &pb_emlrtBCI, (emlrtConstCTX)sp);
              }
              if (((int32_T)(y + 1.0) < 1) ||
                  ((int32_T)(y + 1.0) > numAsters_tmp)) {
                emlrtDynamicBoundsCheckR2012b((int32_T)(y + 1.0), 1,
                                              numAsters_tmp, &ob_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              if (((int32_T)filSubIdx < 1) ||
                  ((int32_T)filSubIdx > astralNum_tmp)) {
                emlrtDynamicBoundsCheckR2012b((int32_T)filSubIdx, 1,
                                              astralNum_tmp, &nb_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              if (((int32_T)nodeCount < 1) ||
                  ((int32_T)nodeCount > nodes->size[0])) {
                emlrtDynamicBoundsCheckR2012b((int32_T)nodeCount, 1,
                                              nodes->size[0], &mb_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              nodes_data[(int32_T)nodeCount - 1] =
                  centers_tmp + l * filSubJdx * t1;
              nodes_data[((int32_T)nodeCount + nodes->size[0]) - 1] =
                  b_centers_tmp + l * denom_tmp * t1;
              if (((int32_T)nodeCount < 1) ||
                  ((int32_T)nodeCount > filCross->size[0])) {
                emlrtDynamicBoundsCheckR2012b((int32_T)nodeCount, 1,
                                              filCross->size[0], &lb_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              filCross_data[(int32_T)nodeCount - 1] = (real_T)idx + 1.0;
              filCross_data[((int32_T)nodeCount + filCross->size[0]) - 1] =
                  b_jdx;
            }
          }
        }
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b((emlrtConstCTX)sp);
        }
      }
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b((emlrtConstCTX)sp);
      }
    }
    emxFree_boolean_T(sp, &tooFar);
  }
  emxFree_real_T(sp, &r1);
  emxFree_real_T(sp, &r);
  emxFree_real_T(sp, &otherAsters);
  emxFree_real_T(sp, &oldNodes);
  emxFree_real_T(sp, &centerSepSQR);
  emxFree_real_T(sp, &sines);
  emxFree_real_T(sp, &cosines);
  if ((int32_T)nodeCount == 0) {
    nodes->size[0] = 0;
    nodes->size[1] = 2;
    filCross->size[0] = 0;
    filCross->size[1] = 2;
  } else {
    if (nodes->size[0] < 1) {
      emlrtDynamicBoundsCheckR2012b(1, 1, nodes->size[0], &kb_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    loop_ub = (int32_T)nodeCount;
    if (((int32_T)nodeCount < 1) || ((int32_T)nodeCount > nodes->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)nodeCount, 1, nodes->size[0],
                                    &jb_emlrtBCI, (emlrtConstCTX)sp);
    }
    for (i = 0; i < 2; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        nodes_data[i1 + (int32_T)nodeCount * i] =
            nodes_data[i1 + nodes->size[0] * i];
      }
    }
    i = nodes->size[0] * nodes->size[1];
    nodes->size[0] = (int32_T)nodeCount;
    nodes->size[1] = 2;
    emxEnsureCapacity_real_T(sp, nodes, i, &hb_emlrtRTEI);
    if (filCross->size[0] < 1) {
      emlrtDynamicBoundsCheckR2012b(1, 1, filCross->size[0], &ib_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)nodeCount < 1) || ((int32_T)nodeCount > filCross->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)nodeCount, 1, filCross->size[0],
                                    &hb_emlrtBCI, (emlrtConstCTX)sp);
    }
    for (i = 0; i < 2; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        filCross_data[i1 + (int32_T)nodeCount * i] =
            filCross_data[i1 + filCross->size[0] * i];
      }
    }
    i = filCross->size[0] * filCross->size[1];
    filCross->size[0] = (int32_T)nodeCount;
    filCross->size[1] = 2;
    emxEnsureCapacity_real_T(sp, filCross, i, &kb_emlrtRTEI);
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void sortNodes(const emlrtStack *sp, real_T filIdx,
                      const emxArray_real_T *nodes,
                      const emxArray_real_T *filCross,
                      const emxArray_real_T *centers, real_T astralNum,
                      emxArray_real_T *nodeOrdering)
{
  __m128d r;
  __m128d r1;
  jmp_buf emlrtJBEnviron;
  jmp_buf *volatile emlrtJBStack;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  emxArray_boolean_T *x;
  emxArray_int32_T *i_tmp;
  emxArray_int32_T *nodesOnFil;
  emxArray_real_T *a;
  emxArray_real_T *dotProducts;
  emxArray_real_T *nodeCoordsOnFil;
  real_T a_tmp[2];
  real_T b_nodesOnFil[2];
  real_T y[2];
  const real_T *centers_data;
  const real_T *filCross_data;
  const real_T *nodes_data;
  real_T d;
  real_T d1;
  real_T d2;
  real_T y_tmp;
  real_T *dotProducts_data;
  real_T *nodeCoordsOnFil_data;
  real_T *nodeOrdering_data;
  int32_T iv[2];
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T idx;
  int32_T ii;
  int32_T jj;
  int32_T loop_ub;
  int32_T loop_ub_tmp;
  int32_T sortNodes_numThreads;
  int32_T *i_tmp_data;
  int32_T *nodesOnFil_data;
  boolean_T emlrtHadParallelError = false;
  boolean_T exitg1;
  boolean_T guard1;
  boolean_T p;
  boolean_T *x_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  centers_data = centers->data;
  filCross_data = filCross->data;
  nodes_data = nodes->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /*  SORTNODES orders the node indices along a filament */
  /*    Inputs: */
  /*        filIdx (scalar): filament index (see note in findNodes for details)
   */
  /*        nodes (numNodes x 2 double): list of (x,y) coordinates of filament
   */
  /*        crossings (EXCLUDING astral centers) */
  /*        filCross (numNodes x 2 double): list of pairs of filament indices */
  /*        corresponding to the filaments that cross at a particular node; */
  /*        listed so that filCross(idx,1) < filCross(idx,2) */
  /*        centers (numAsters x 2 double): (x,y) coordinates of astral centers
   */
  /*        astralNum (scalar): number of filaments per aster */
  /*    Outputs: */
  /*        nodeOrdering (double vector): sorted list of node indices appearing
   */
  /*        on the specified filament. First element is the node closest to the
   */
  /*        astral center. If there are no nodes on this filament, returns [].
   */
  st.site = &tb_emlrtRSI;
  emxInit_boolean_T(&st, &x, 2, &fc_emlrtRTEI);
  i = x->size[0] * x->size[1];
  x->size[0] = filCross->size[0];
  x->size[1] = 2;
  emxEnsureCapacity_boolean_T(&st, x, i, &fc_emlrtRTEI);
  x_data = x->data;
  loop_ub_tmp = filCross->size[0] << 1;
  if (loop_ub_tmp < 3200) {
    for (i1 = 0; i1 < loop_ub_tmp; i1++) {
      x_data[i1] = (filCross_data[i1] == filIdx);
    }
  } else {
    emlrtEnterParallelRegion(&st, omp_in_parallel());
    emlrtPushJmpBuf(&st, &emlrtJBStack);
    sortNodes_numThreads = emlrtAllocRegionTLSs(
        st.tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(sortNodes_numThreads)

    for (i1 = 0; i1 < loop_ub_tmp; i1++) {
      x_data[i1] = (filCross_data[i1] == filIdx);
    }
    emlrtPopJmpBuf(&st, &emlrtJBStack);
    emlrtExitParallelRegion(&st, omp_in_parallel());
  }
  b_st.site = &yb_emlrtRSI;
  emxInit_int32_T(&b_st, &i_tmp, &rc_emlrtRTEI);
  i_tmp_data = i_tmp->data;
  if (loop_ub_tmp == 0) {
    i_tmp->size[0] = 0;
  } else {
    c_st.site = &bc_emlrtRSI;
    idx = 0;
    i = i_tmp->size[0];
    i_tmp->size[0] = loop_ub_tmp;
    emxEnsureCapacity_int32_T(&c_st, i_tmp, i, &gc_emlrtRTEI);
    i_tmp_data = i_tmp->data;
    ii = 1;
    jj = 1;
    exitg1 = false;
    while ((!exitg1) && (jj <= 2)) {
      guard1 = false;
      if (x_data[(ii + x->size[0] * (jj - 1)) - 1]) {
        idx++;
        i_tmp_data[idx - 1] = ii;
        if (idx >= loop_ub_tmp) {
          exitg1 = true;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }
      if (guard1) {
        ii++;
        if (ii > x->size[0]) {
          ii = 1;
          jj++;
        }
      }
    }
    if (idx > loop_ub_tmp) {
      emlrtErrorWithMessageIdR2018a(&c_st, &k_emlrtRTEI,
                                    "Coder:builtins:AssertionFailed",
                                    "Coder:builtins:AssertionFailed", 0);
    }
    if (loop_ub_tmp == 1) {
      if (idx == 0) {
        i_tmp->size[0] = 0;
      }
    } else {
      if (idx < 1) {
        i = 0;
      } else {
        i = idx;
      }
      iv[0] = 1;
      iv[1] = i;
      d_st.site = &cc_emlrtRSI;
      indexShapeCheck(&d_st, i_tmp->size[0], iv);
      ii = i_tmp->size[0];
      i_tmp->size[0] = i;
      emxEnsureCapacity_int32_T(&c_st, i_tmp, ii, &ic_emlrtRTEI);
      i_tmp_data = i_tmp->data;
      iv[0] = 1;
      iv[1] = i;
      d_st.site = &dc_emlrtRSI;
      indexShapeCheck(&d_st, loop_ub_tmp, iv);
      iv[0] = 1;
      iv[1] = i;
      d_st.site = &ec_emlrtRSI;
      indexShapeCheck(&d_st, loop_ub_tmp, iv);
    }
  }
  if ((x->size[0] == 1) && (i_tmp->size[0] != 1)) {
    c_st.site = &ac_emlrtRSI;
    warning(&c_st);
  }
  emxFree_boolean_T(&b_st, &x);
  emxInit_int32_T(&st, &nodesOnFil, &qc_emlrtRTEI);
  loop_ub = i_tmp->size[0];
  i = nodesOnFil->size[0];
  nodesOnFil->size[0] = i_tmp->size[0];
  emxEnsureCapacity_int32_T(&st, nodesOnFil, i, &hc_emlrtRTEI);
  nodesOnFil_data = nodesOnFil->data;
  loop_ub_tmp = i_tmp->size[0];
  if (i_tmp->size[0] < 3200) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      nodesOnFil_data[i1] = i_tmp_data[i1];
    }
  } else {
    emlrtEnterParallelRegion(&st, omp_in_parallel());
    emlrtPushJmpBuf(&st, &emlrtJBStack);
    sortNodes_numThreads = emlrtAllocRegionTLSs(
        st.tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(sortNodes_numThreads)

    for (i1 = 0; i1 < loop_ub_tmp; i1++) {
      nodesOnFil_data[i1] = i_tmp_data[i1];
    }
    emlrtPopJmpBuf(&st, &emlrtJBStack);
    emlrtExitParallelRegion(&st, omp_in_parallel());
  }
  emxInit_real_T(sp, &nodeCoordsOnFil, 2, &kc_emlrtRTEI);
  emxInit_real_T(sp, &dotProducts, 1, &pc_emlrtRTEI);
  emxInit_real_T(sp, &a, 2, &sc_emlrtRTEI);
  if (nodesOnFil->size[0] == 0) {
    nodeOrdering->size[0] = 0;
    nodeOrdering->size[1] = 0;
  } else if (nodesOnFil->size[0] == 1) {
    i = nodeOrdering->size[0] * nodeOrdering->size[1];
    nodeOrdering->size[0] = 1;
    nodeOrdering->size[1] = 1;
    emxEnsureCapacity_real_T(sp, nodeOrdering, i, &jc_emlrtRTEI);
    nodeOrdering_data = nodeOrdering->data;
    nodeOrdering_data[0] = nodesOnFil_data[0];
  } else {
    ii = nodes->size[0];
    i = nodeCoordsOnFil->size[0] * nodeCoordsOnFil->size[1];
    nodeCoordsOnFil->size[0] = loop_ub;
    nodeCoordsOnFil->size[1] = 2;
    emxEnsureCapacity_real_T(sp, nodeCoordsOnFil, i, &kc_emlrtRTEI);
    nodeCoordsOnFil_data = nodeCoordsOnFil->data;
    if ((nodesOnFil->size[0] << 1) < 3200) {
      for (i1 = 0; i1 < 2; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          if ((nodesOnFil_data[i2] < 1) || (nodesOnFil_data[i2] > ii)) {
            emlrtDynamicBoundsCheckR2012b(nodesOnFil_data[i2], 1, ii,
                                          &rd_emlrtBCI, (emlrtConstCTX)sp);
          }
          nodeCoordsOnFil_data[i2 + nodeCoordsOnFil->size[0] * i1] =
              nodes_data[(nodesOnFil_data[i2] + nodes->size[0] * i1) - 1];
        }
      }
    } else {
      emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
      emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      sortNodes_numThreads =
          emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(),
                               omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel num_threads(sortNodes_numThreads) private(                \
        i2, e_st, emlrtJBEnviron) firstprivate(emlrtHadParallelError)
      {
        if (setjmp(emlrtJBEnviron) == 0) {
          e_st.prev = sp;
          e_st.tls = emlrtAllocTLS((emlrtCTX)sp, omp_get_thread_num());
          e_st.site = NULL;
          emlrtSetJmpBuf(&e_st, &emlrtJBEnviron);
        } else {
          emlrtHadParallelError = true;
        }
#pragma omp for nowait
        for (i1 = 0; i1 < 2; i1++) {
          if (emlrtHadParallelError) {
            continue;
          }
          if (setjmp(emlrtJBEnviron) == 0) {
            for (i2 = 0; i2 < loop_ub; i2++) {
              if ((nodesOnFil_data[i2] < 1) || (nodesOnFil_data[i2] > ii)) {
                emlrtDynamicBoundsCheckR2012b(nodesOnFil_data[i2], 1, ii,
                                              &rd_emlrtBCI, &e_st);
              }
              nodeCoordsOnFil_data[i2 + nodeCoordsOnFil->size[0] * i1] =
                  nodes_data[(nodesOnFil_data[i2] + nodes->size[0] * i1) - 1];
            }
          } else {
            emlrtHadParallelError = true;
          }
        }
      }
      emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
    }
    y_tmp = muDoubleScalarFloor((filIdx - 1.0) / astralNum);
    if (((int32_T)(y_tmp + 1.0) < 1) ||
        ((int32_T)(y_tmp + 1.0) > centers->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(y_tmp + 1.0), 1, centers->size[0],
                                    &qd_emlrtBCI, (emlrtConstCTX)sp);
    }
    d = centers_data[(int32_T)(y_tmp + 1.0) - 1];
    d1 = nodeCoordsOnFil_data[0] - d;
    a_tmp[0] = d1;
    y[0] = d1 * d1;
    d2 = centers_data[((int32_T)(y_tmp + 1.0) + centers->size[0]) - 1];
    d1 = nodeCoordsOnFil_data[nodeCoordsOnFil->size[0]] - d2;
    a_tmp[1] = d1;
    y[1] = d1 * d1;
    y_tmp = sumColumnB(y);
    st.site = &ub_emlrtRSI;
    if (y_tmp < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &l_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    y_tmp = muDoubleScalarSqrt(y_tmp);
    /*  next two lines find dot product "direction * (node - r0)" for all nodes
     */
    y[0] = d;
    y[1] = d2;
    idx = nodesOnFil->size[0];
    b_nodesOnFil[0] = nodesOnFil->size[0];
    b_nodesOnFil[1] = 1.0;
    st.site = &vb_emlrtRSI;
    repmat(&st, y, b_nodesOnFil, a);
    nodeOrdering_data = a->data;
    if ((loop_ub != a->size[0]) && (a->size[0] != 1)) {
      emlrtDimSizeImpxCheckR2021b(loop_ub, a->size[0], &h_emlrtECI,
                                  (emlrtConstCTX)sp);
    }
    loop_ub_tmp = nodeCoordsOnFil->size[0] << 1;
    i = nodeCoordsOnFil->size[0] * nodeCoordsOnFil->size[1];
    nodeCoordsOnFil->size[1] = 2;
    emxEnsureCapacity_real_T(sp, nodeCoordsOnFil, i, &lc_emlrtRTEI);
    nodeCoordsOnFil_data = nodeCoordsOnFil->data;
    ii = (loop_ub_tmp / 2) << 1;
    jj = ii - 2;
    for (i1 = 0; i1 <= jj; i1 += 2) {
      r = _mm_loadu_pd(&nodeCoordsOnFil_data[i1]);
      r1 = _mm_loadu_pd(&nodeOrdering_data[i1]);
      _mm_storeu_pd(&nodeCoordsOnFil_data[i1], _mm_sub_pd(r, r1));
    }
    if (loop_ub_tmp - ii < 3200) {
      for (i1 = ii; i1 < loop_ub_tmp; i1++) {
        nodeCoordsOnFil_data[i1] -= nodeOrdering_data[i1];
      }
    } else {
      emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
      emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      sortNodes_numThreads =
          emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(),
                               omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(sortNodes_numThreads)

      for (i1 = ii; i1 < loop_ub_tmp; i1++) {
        nodeCoordsOnFil_data[i1] -= nodeOrdering_data[i1];
      }
      emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
    }
    st.site = &wb_emlrtRSI;
    r = _mm_loadu_pd(&a_tmp[0]);
    _mm_storeu_pd(&a_tmp[0], _mm_div_pd(r, _mm_set1_pd(y_tmp)));
    b_nodesOnFil[0] = idx;
    b_nodesOnFil[1] = 1.0;
    b_st.site = &wb_emlrtRSI;
    repmat(&b_st, a_tmp, b_nodesOnFil, a);
    nodeOrdering_data = a->data;
    a_tmp[0] = a->size[0];
    y[0] = nodeCoordsOnFil->size[0];
    a_tmp[1] = 2.0;
    y[1] = 2.0;
    p = true;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii < 2)) {
      if (a_tmp[ii] != y[ii]) {
        p = false;
        exitg1 = true;
      } else {
        ii++;
      }
    }
    if (!p) {
      emlrtErrorWithMessageIdR2018a(&st, &m_emlrtRTEI,
                                    "MATLAB:dot:InputSizeMismatch",
                                    "MATLAB:dot:InputSizeMismatch", 0);
    }
    i = a->size[0];
    ii = dotProducts->size[0];
    dotProducts->size[0] = a->size[0];
    emxEnsureCapacity_real_T(&st, dotProducts, ii, &mc_emlrtRTEI);
    dotProducts_data = dotProducts->data;
    ii = -1;
    b_st.site = &gc_emlrtRSI;
    if (a->size[0] > 2147483646) {
      c_st.site = &t_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (idx = 0; idx < i; idx++) {
      ii++;
      jj = ii + i;
      dotProducts_data[idx] = nodeOrdering_data[ii] * nodeCoordsOnFil_data[ii] +
                              nodeOrdering_data[jj] * nodeCoordsOnFil_data[jj];
    }
    /*  now sort the dot products in increasing order and sort the node indices
     */
    /*  using the same permutation */
    st.site = &xb_emlrtRSI;
    b_st.site = &hc_emlrtRSI;
    sort(&b_st, dotProducts, i_tmp);
    i_tmp_data = i_tmp->data;
    loop_ub_tmp = i_tmp->size[0];
    i = dotProducts->size[0];
    dotProducts->size[0] = i_tmp->size[0];
    emxEnsureCapacity_real_T(sp, dotProducts, i, &nc_emlrtRTEI);
    dotProducts_data = dotProducts->data;
    ii = i_tmp->size[0];
    i = (i_tmp->size[0] < 3200);
    if (i) {
      for (i1 = 0; i1 < loop_ub_tmp; i1++) {
        if ((i_tmp_data[i1] < 1) || (i_tmp_data[i1] > loop_ub)) {
          emlrtDynamicBoundsCheckR2012b(i_tmp_data[i1], 1, loop_ub,
                                        &sd_emlrtBCI, (emlrtConstCTX)sp);
        }
        dotProducts_data[i1] = nodesOnFil_data[i_tmp_data[i1] - 1];
      }
    } else {
      emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
      emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      sortNodes_numThreads =
          emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(),
                               omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel num_threads(sortNodes_numThreads) private(                \
        emlrtJBEnviron, f_st) firstprivate(emlrtHadParallelError)
      {
        if (setjmp(emlrtJBEnviron) == 0) {
          f_st.prev = sp;
          f_st.tls = emlrtAllocTLS((emlrtCTX)sp, omp_get_thread_num());
          f_st.site = NULL;
          emlrtSetJmpBuf(&f_st, &emlrtJBEnviron);
        } else {
          emlrtHadParallelError = true;
        }
#pragma omp for nowait
        for (i1 = 0; i1 < ii; i1++) {
          if (emlrtHadParallelError) {
            continue;
          }
          if (setjmp(emlrtJBEnviron) == 0) {
            if ((i_tmp_data[i1] < 1) || (i_tmp_data[i1] > loop_ub)) {
              emlrtDynamicBoundsCheckR2012b(i_tmp_data[i1], 1, loop_ub,
                                            &sd_emlrtBCI, &f_st);
            }
            dotProducts_data[i1] = nodesOnFil_data[i_tmp_data[i1] - 1];
          } else {
            emlrtHadParallelError = true;
          }
        }
      }
      emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
    }
    ii = nodeOrdering->size[0] * nodeOrdering->size[1];
    nodeOrdering->size[0] = loop_ub_tmp;
    nodeOrdering->size[1] = 1;
    emxEnsureCapacity_real_T(sp, nodeOrdering, ii, &oc_emlrtRTEI);
    nodeOrdering_data = nodeOrdering->data;
    if (i) {
      for (i1 = 0; i1 < loop_ub_tmp; i1++) {
        nodeOrdering_data[i1] = dotProducts_data[i1];
      }
    } else {
      emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
      emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      sortNodes_numThreads =
          emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(),
                               omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(sortNodes_numThreads)

      for (i1 = 0; i1 < loop_ub_tmp; i1++) {
        nodeOrdering_data[i1] = dotProducts_data[i1];
      }
      emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
    }
  }
  emxFree_real_T(sp, &a);
  emxFree_int32_T(sp, &i_tmp);
  emxFree_int32_T(sp, &nodesOnFil);
  emxFree_real_T(sp, &dotProducts);
  emxFree_real_T(sp, &nodeCoordsOnFil);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

emlrtCTX emlrtGetRootTLSGlobal(void)
{
  return emlrtRootTLSGlobal;
}

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, emlrtConstCTX aTLS,
                         void *aData)
{
  omp_set_lock(&emlrtLockGlobal);
  emlrtCallLockeeFunction(aLockee, aTLS, aData);
  omp_unset_lock(&emlrtLockGlobal);
}

void generateAstralNetwork(const emlrtStack *sp, real_T numAsters, real_T l,
                           real_T D, real_T astralNum, boolean_T nodesOnly,
                           struct0_T *network, struct1_T *crossings,
                           struct2_T *asters)
{
  __m128d r;
  jmp_buf *volatile emlrtJBStack;
  emlrtStack b_st;
  emlrtStack st;
  emxArray_real_T *b;
  emxArray_real_T *r1;
  emxArray_real_T *springs;
  emxArray_real_T *thisOrder;
  real_T b_numAsters[2];
  real_T numFil;
  real_T *springs_data;
  real_T *thisOrder_data;
  int32_T a[2];
  int32_T newSpring_size[2];
  int32_T generateAstralNetwork_numThreads;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T jdx;
  int32_T loop_ub;
  int32_T maxdimlen;
  int32_T n;
  int32_T vectorUB;
  boolean_T out;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /*  GENERATEASTRALNETWORK Constructs an astral network and reports properties
   */
  /*    Inputs: */
  /*        numAsters (scalar): (whole) number of asters to distribute */
  /*        l (scalar): length of individual filament */
  /*        D (scalar): domain size; asters are distributed in the square */
  /*        with corners at (0,0) and (D,D) */
  /*        astralNum (scalar): (whole) number of filaments per aster */
  /*        nodesOnly (boolean): passing true generates nodes only and skips */
  /*        creation of network fields augNodes, springs, & ends, pass false to
   */
  /*        generate all fields */
  /*    Outputs: */
  /*        network (struct): has fields 'nodes', 'augNodes', 'springs', 'ends'
   */
  /*            'nodes': see findNodes auxiliary function */
  /*            'augNodes': see defineSprings auxiliary function */
  /*            'springs': see defineSprings auxiliary function */
  /*            'ends': see defineSprings auxiliary function */
  /*        crossings (struct): has fields 'filCross' and 'centerCross' */
  /*            'filCross': see findNodes auxiliary function */
  /*            'centerCross': which filaments cross at each astral center; */
  /*            each row represents a center */
  /*        asters (struct): has fields 'centers', 'orients' */
  /*            'centers': (x,y) coordinates of astral centers */
  /*            'orients': angles giving each filament's orientation about its
   */
  /*            astral center */
  /*    Note: currently no boundary/edge is constructed, and the implementation
   */
  /*    does not store/compile information for running Metropolis-Hastings */
  if (numAsters == 0.0) {
    st.site = &emlrtRSI;
    b_error(&st);
  }
  b_numAsters[0] = numAsters;
  b_numAsters[1] = 2.0;
  st.site = &b_emlrtRSI;
  b_rand(&st, b_numAsters, network->nodes);
  i = asters->centers->size[0] * asters->centers->size[1];
  asters->centers->size[0] = network->nodes->size[0];
  asters->centers->size[1] = 2;
  emxEnsureCapacity_real_T(sp, asters->centers, i, &p_emlrtRTEI);
  loop_ub = network->nodes->size[0] << 1;
  maxdimlen = (loop_ub / 2) << 1;
  vectorUB = maxdimlen - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&network->nodes->data[i]);
    _mm_storeu_pd(&asters->centers->data[i], _mm_mul_pd(_mm_set1_pd(D), r));
  }
  for (i = maxdimlen; i < loop_ub; i++) {
    asters->centers->data[i] = D * network->nodes->data[i];
  }
  b_numAsters[0] = numAsters;
  b_numAsters[1] = astralNum;
  emxInit_real_T(sp, &b, 2, &ab_emlrtRTEI);
  st.site = &c_emlrtRSI;
  c_rand(&st, b_numAsters, b);
  thisOrder_data = b->data;
  i = asters->orients->size[0] * asters->orients->size[1];
  asters->orients->size[0] = b->size[0];
  asters->orients->size[1] = b->size[1];
  emxEnsureCapacity_real_T(sp, asters->orients, i, &q_emlrtRTEI);
  loop_ub = b->size[0] * b->size[1];
  maxdimlen = (loop_ub / 2) << 1;
  vectorUB = maxdimlen - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&thisOrder_data[i]);
    _mm_storeu_pd(&asters->orients->data[i],
                  _mm_mul_pd(_mm_set1_pd(6.2831853071795862), r));
  }
  for (i = maxdimlen; i < loop_ub; i++) {
    asters->orients->data[i] = 6.2831853071795862 * thisOrder_data[i];
  }
  emxFree_real_T(sp, &b);
  st.site = &d_emlrtRSI;
  findNodes(&st, asters->centers, asters->orients, l, network->nodes,
            crossings->filCross);
  numFil = numAsters * astralNum;
  /*  aster idx 1 groups filaments 1,2,...,astralNum; and so on */
  emxInit_real_T(sp, &r1, 2, &y_emlrtRTEI);
  thisOrder_data = r1->data;
  if (muDoubleScalarIsNaN(numFil)) {
    i = r1->size[0] * r1->size[1];
    r1->size[0] = 1;
    r1->size[1] = 1;
    emxEnsureCapacity_real_T(sp, r1, i, &r_emlrtRTEI);
    thisOrder_data = r1->data;
    thisOrder_data[0] = rtNaN;
  } else if (numFil < 1.0) {
    r1->size[0] = 1;
    r1->size[1] = 0;
  } else {
    i = r1->size[0] * r1->size[1];
    r1->size[0] = 1;
    r1->size[1] = (int32_T)(numFil - 1.0) + 1;
    emxEnsureCapacity_real_T(sp, r1, i, &r_emlrtRTEI);
    thisOrder_data = r1->data;
    loop_ub = (int32_T)(numFil - 1.0);
    if ((int32_T)(numFil - 1.0) + 1 < 3200) {
      for (i1 = 0; i1 <= loop_ub; i1++) {
        thisOrder_data[i1] = (real_T)i1 + 1.0;
      }
    } else {
      emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
      emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      generateAstralNetwork_numThreads =
          emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(),
                               omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(generateAstralNetwork_numThreads)

      for (i1 = 0; i1 <= loop_ub; i1++) {
        thisOrder_data[i1] = (real_T)i1 + 1.0;
      }
      emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
      emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
    }
  }
  b_numAsters[0] = astralNum;
  b_numAsters[1] = numAsters;
  st.site = &e_emlrtRSI;
  maxdimlen = r1->size[1];
  b_st.site = &jb_emlrtRSI;
  assertValidSizeArg(&b_st, b_numAsters);
  vectorUB = (int32_T)astralNum;
  loop_ub = (int32_T)numAsters;
  n = 1;
  if (r1->size[1] > 1) {
    n = r1->size[1];
  }
  maxdimlen = muIntScalarMax_sint32(maxdimlen, n);
  if ((int32_T)astralNum > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if ((int32_T)numAsters > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  out = ((int32_T)astralNum >= 0);
  if ((!out) || ((int32_T)numAsters < 0)) {
    out = false;
  }
  if (!out) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "MATLAB:checkDimCommon:nonnegativeSize",
                                  "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }
  if ((int32_T)astralNum * (int32_T)numAsters != r1->size[1]) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  i = crossings->centerCross->size[0] * crossings->centerCross->size[1];
  crossings->centerCross->size[0] = (int32_T)numAsters;
  crossings->centerCross->size[1] = (int32_T)astralNum;
  emxEnsureCapacity_real_T(sp, crossings->centerCross, i, &s_emlrtRTEI);
  for (i = 0; i < vectorUB; i++) {
    for (i2 = 0; i2 < loop_ub; i2++) {
      crossings->centerCross->data[i2 + crossings->centerCross->size[0] * i] =
          thisOrder_data[i + (int32_T)astralNum * i2];
    }
  }
  emxFree_real_T(sp, &r1);
  if (nodesOnly) {
    /*  skip spring definition */
    network->augNodes->size[0] = 0;
    network->augNodes->size[1] = 2;
    network->springs->size[0] = 0;
    network->springs->size[1] = 4;
    network->ends->size[0] = 0;
    network->ends->size[1] = 2;
  } else {
    real_T numSpringsGuess;
    int32_T c_numAsters;
    int32_T numAsters_tmp;
    uint32_T springCount;
    /*  full network generation */
    st.site = &f_emlrtRSI;
    /*  DEFINESPRINGS partitions filaments into node-bounded sub-segments  */
    /*  (springs) and dangling ends */
    /*    Inputs: */
    /*        nodes (numNodes x 2 double): list of (x,y) coordinates of filament
     */
    /*        crossings (EXCLUDING astral centers) */
    /*        filCross (numNodes x 2 double): list of pairs of filament indices
     */
    /*        corresponding to the filaments that cross at a particular node; */
    /*        listed so that filCross(idx,1) < filCross(idx,2) */
    /*        centers (numAsters x 2 double): (x,y) coordinates of astral
     * centers */
    /*        l (scalar): length of individual filament */
    /*        astralNum (scalar): number of filaments per aster */
    /*    Outputs: */
    /*        augNodes (numAsters+numNodes x 2 double)**: conditionally
     * augmented  */
    /*        node list; if astralNum >= 2, rows 1:numAsters contain (x,y)  */
    /*        coordinates of astral centers, remaining rows contain (x,y)  */
    /*        coordinates of inter-aster nodes */
    /*        **NOTE: if astralNum == 1, augNodes is the same as nodes (centers
     */
    /*        are fictitious), so its size is (numNodes x 2) */
    /*        springs (numSprings x 4 double): descriptions of segments bounded
     */
    /*        by nodes at each end, each row has the structure */
    /*            (1) lesser augNode index at one end of the spring */
    /*            (2) greater augNode index at the other end of the spring */
    /*            (3) filament index (which filament is the spring on) */
    /*            (4) the original distance between the nodes */
    /*        ends (numFil x 2 double): lengths of segments on the ends of each
     */
    /*        filament. ends(:,2) lists the "right" dangling ends, distal to the
     */
    /*        astral center. ends(:,1) lists the "left" dangling ends, proximal
     */
    /*        to the astral center; these are all 0 if astralNum >= 2. */
    numAsters_tmp = asters->centers->size[0];
    c_numAsters = asters->centers->size[0];
    numFil = (real_T)asters->centers->size[0] * astralNum;
    numSpringsGuess = 10.0 * (real_T)network->nodes->size[0];
    emxInit_real_T(&st, &springs, 2, &x_emlrtRTEI);
    if (numSpringsGuess != (int32_T)numSpringsGuess) {
      emlrtIntegerCheckR2012b(numSpringsGuess, &e_emlrtDCI, &st);
    }
    i = springs->size[0] * springs->size[1];
    springs->size[0] = (int32_T)numSpringsGuess;
    springs->size[1] = 4;
    emxEnsureCapacity_real_T(&st, springs, i, &t_emlrtRTEI);
    springs_data = springs->data;
    if (numSpringsGuess != (int32_T)numSpringsGuess) {
      emlrtIntegerCheckR2012b(numSpringsGuess, &b_emlrtDCI, &st);
    }
    loop_ub = (int32_T)numSpringsGuess << 2;
    if (loop_ub < 3200) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        springs_data[i1] = 0.0;
      }
    } else {
      emlrtEnterParallelRegion(&st, omp_in_parallel());
      emlrtPushJmpBuf(&st, &emlrtJBStack);
      generateAstralNetwork_numThreads =
          emlrtAllocRegionTLSs(st.tls, omp_in_parallel(), omp_get_max_threads(),
                               omp_get_num_procs());
#pragma omp parallel for num_threads(generateAstralNetwork_numThreads)

      for (i1 = 0; i1 < loop_ub; i1++) {
        springs_data[i1] = 0.0;
      }
      emlrtPopJmpBuf(&st, &emlrtJBStack);
      emlrtExitParallelRegion(&st, omp_in_parallel());
    }
    if (!(numFil >= 0.0)) {
      emlrtNonNegativeCheckR2012b(numFil, &c_emlrtDCI, &st);
    }
    i = (int32_T)muDoubleScalarFloor(numFil);
    if (numFil != i) {
      emlrtIntegerCheckR2012b(numFil, &d_emlrtDCI, &st);
    }
    i2 = (int32_T)numFil;
    i3 = network->ends->size[0] * network->ends->size[1];
    network->ends->size[0] = i2;
    network->ends->size[1] = 2;
    emxEnsureCapacity_real_T(&st, network->ends, i3, &u_emlrtRTEI);
    if (i2 != i) {
      emlrtIntegerCheckR2012b(numFil, &emlrtDCI, &st);
    }
    loop_ub = i2 << 1;
    for (i = 0; i < loop_ub; i++) {
      network->ends->data[i] = 0.0;
    }
    springCount = 0U;
    network->augNodes->size[0] = 0;
    network->augNodes->size[1] = 2;
    emxInit_real_T(&st, &thisOrder, 2, &bb_emlrtRTEI);
    if (astralNum == 1.0) {
      /*  routine for "Classical Mikado" networks */
      i = network->nodes->size[0];
      i3 = network->augNodes->size[0] * network->augNodes->size[1];
      network->augNodes->size[0] = network->nodes->size[0];
      network->augNodes->size[1] = 2;
      emxEnsureCapacity_real_T(&st, network->augNodes, i3, &v_emlrtRTEI);
      loop_ub = network->nodes->size[0] << 1;
      for (i3 = 0; i3 < loop_ub; i3++) {
        network->augNodes->data[i3] = network->nodes->data[i3];
      }
      /*  centers are fictitious */
      emlrtForLoopVectorCheckR2021a(1.0, 1.0, numFil, mxDOUBLE_CLASS,
                                    (int32_T)numFil, &e_emlrtRTEI, &st);
      for (loop_ub = 0; loop_ub < i2; loop_ub++) {
        b_st.site = &sb_emlrtRSI;
        sortNodes(&b_st, (real_T)loop_ub + 1.0, network->nodes,
                  crossings->filCross, asters->centers, 1.0, thisOrder);
        thisOrder_data = thisOrder->data;
        if ((thisOrder->size[0] == 0) || (thisOrder->size[1] == 0)) {
          /*  0 nodes = filament doesn't touch any other asters! */
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &i_emlrtBCI,
                                          &st);
          }
          network->ends->data[loop_ub + network->ends->size[0]] = l;
          /*  treat unused filament as "right" dangling end */
        } else if (thisOrder->size[0] == 1) {
          /*  1 node + astralNum == 1 -> two dangling ends */
          if (loop_ub + 1 > numAsters_tmp) {
            emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1, numAsters_tmp,
                                          &db_emlrtBCI, &st);
          }
          if (((int32_T)thisOrder_data[0] < 1) ||
              ((int32_T)thisOrder_data[0] > i)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)thisOrder_data[0], 1, i,
                                          &cb_emlrtBCI, &st);
          }
          numFil = network->nodes->data[(int32_T)thisOrder_data[0] - 1] -
                   asters->centers->data[loop_ub];
          b_numAsters[0] = numFil * numFil;
          numFil = network->nodes->data[((int32_T)thisOrder_data[0] +
                                         network->nodes->size[0]) -
                                        1] -
                   asters->centers->data[loop_ub + asters->centers->size[0]];
          b_numAsters[1] = numFil * numFil;
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &h_emlrtBCI,
                                          &st);
          }
          network->ends->data[loop_ub] = sumColumnB(b_numAsters);
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &h_emlrtBCI,
                                          &st);
          }
          b_st.site = &rb_emlrtRSI;
          thisOrder_data = &network->ends->data[loop_ub];
          if (*thisOrder_data < 0.0) {
            emlrtErrorWithMessageIdR2018a(
                &b_st, &l_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
          }
          *thisOrder_data = muDoubleScalarSqrt(*thisOrder_data);
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &f_emlrtBCI,
                                          &st);
          }
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &g_emlrtBCI,
                                          &st);
          }
          network->ends->data[loop_ub + network->ends->size[0]] =
              l - network->ends->data[loop_ub + network->ends->size[0]];
        } else {
          real_T totSpringLength;
          /*  >=2 nodes = there are certainly springs along this filament */
          if (loop_ub + 1 > numAsters_tmp) {
            emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1, numAsters_tmp,
                                          &bb_emlrtBCI, &st);
          }
          totSpringLength = 0.0;
          /*  pair adjacent nodes to form springs */
          n = thisOrder->size[0] - 2;
          maxdimlen = thisOrder->size[0];
          for (jdx = 0; jdx <= n; jdx++) {
            real_T springLength;
            springCount++;
            if (jdx + 1 > maxdimlen) {
              emlrtDynamicBoundsCheckR2012b(jdx + 1, 1, maxdimlen, &ab_emlrtBCI,
                                            &st);
            }
            i3 = (int32_T)thisOrder_data[jdx];
            if ((i3 < 1) || (i3 > i)) {
              emlrtDynamicBoundsCheckR2012b(i3, 1, i, &y_emlrtBCI, &st);
            }
            if (jdx + 2 > maxdimlen) {
              emlrtDynamicBoundsCheckR2012b(jdx + 2, 1, maxdimlen, &x_emlrtBCI,
                                            &st);
            }
            vectorUB = (int32_T)thisOrder_data[jdx + 1];
            if ((vectorUB < 1) || (vectorUB > i)) {
              emlrtDynamicBoundsCheckR2012b(vectorUB, 1, i, &w_emlrtBCI, &st);
            }
            numFil = network->nodes->data[i3 - 1] -
                     network->nodes->data[vectorUB - 1];
            b_numAsters[0] = numFil * numFil;
            numFil =
                network->nodes->data[(i3 + network->nodes->size[0]) - 1] -
                network->nodes->data[(vectorUB + network->nodes->size[0]) - 1];
            b_numAsters[1] = numFil * numFil;
            springLength = sumColumnB(b_numAsters);
            b_st.site = &qb_emlrtRSI;
            if (springLength < 0.0) {
              emlrtErrorWithMessageIdR2018a(
                  &b_st, &l_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                  "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
            }
            springLength = muDoubleScalarSqrt(springLength);
            numSpringsGuess = thisOrder_data[jdx + 1];
            if ((thisOrder_data[jdx] > numSpringsGuess) ||
                (muDoubleScalarIsNaN(thisOrder_data[jdx]) &&
                 (!muDoubleScalarIsNaN(numSpringsGuess)))) {
              numFil = thisOrder_data[jdx + 1];
            } else {
              numFil = thisOrder_data[jdx];
            }
            if ((!(thisOrder_data[jdx] < numSpringsGuess)) &&
                ((!muDoubleScalarIsNaN(thisOrder_data[jdx])) ||
                 muDoubleScalarIsNaN(numSpringsGuess))) {
              numSpringsGuess = thisOrder_data[jdx];
            }
            if (((int32_T)springCount < 1) ||
                ((int32_T)springCount > springs->size[0])) {
              emlrtDynamicBoundsCheckR2012b(
                  (int32_T)springCount, 1, springs->size[0], &fb_emlrtBCI, &st);
            }
            springs_data[(int32_T)springCount - 1] = numFil;
            springs_data[((int32_T)springCount + springs->size[0]) - 1] =
                numSpringsGuess;
            springs_data[((int32_T)springCount + springs->size[0] * 2) - 1] =
                (real_T)loop_ub + 1.0;
            springs_data[((int32_T)springCount + springs->size[0] * 3) - 1] =
                springLength;
            totSpringLength += springLength;
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b(&st);
            }
          }
          if (((int32_T)thisOrder_data[0] < 1) ||
              ((int32_T)thisOrder_data[0] > i)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)thisOrder_data[0], 1, i,
                                          &v_emlrtBCI, &st);
          }
          numFil = network->nodes->data[(int32_T)thisOrder_data[0] - 1] -
                   asters->centers->data[loop_ub];
          b_numAsters[0] = numFil * numFil;
          numFil = network->nodes->data[((int32_T)thisOrder_data[0] +
                                         network->nodes->size[0]) -
                                        1] -
                   asters->centers->data[loop_ub + asters->centers->size[0]];
          b_numAsters[1] = numFil * numFil;
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &e_emlrtBCI,
                                          &st);
          }
          network->ends->data[loop_ub] = sumColumnB(b_numAsters);
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &e_emlrtBCI,
                                          &st);
          }
          b_st.site = &pb_emlrtRSI;
          thisOrder_data = &network->ends->data[loop_ub];
          if (*thisOrder_data < 0.0) {
            emlrtErrorWithMessageIdR2018a(
                &b_st, &l_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
          }
          *thisOrder_data = muDoubleScalarSqrt(*thisOrder_data);
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &c_emlrtBCI,
                                          &st);
          }
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &d_emlrtBCI,
                                          &st);
          }
          network->ends->data[loop_ub + network->ends->size[0]] =
              (l - totSpringLength) - network->ends->data[loop_ub];
        }
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(&st);
        }
      }
    } else if (astralNum >= 2.0) {
      /*  routine for "Astral Mikado" networks */
      i = network->augNodes->size[0] * network->augNodes->size[1];
      network->augNodes->size[0] =
          asters->centers->size[0] + network->nodes->size[0];
      network->augNodes->size[1] = 2;
      emxEnsureCapacity_real_T(&st, network->augNodes, i, &v_emlrtRTEI);
      for (i = 0; i < 2; i++) {
        for (i3 = 0; i3 < numAsters_tmp; i3++) {
          network->augNodes->data[i3 + network->augNodes->size[0] * i] =
              asters->centers->data[i3 + asters->centers->size[0] * i];
        }
      }
      maxdimlen = network->nodes->size[0];
      for (i = 0; i < 2; i++) {
        for (i3 = 0; i3 < maxdimlen; i3++) {
          network->augNodes->data[(i3 + asters->centers->size[0]) +
                                  network->augNodes->size[0] * i] =
              network->nodes->data[i3 + network->nodes->size[0] * i];
        }
      }
      emlrtForLoopVectorCheckR2021a(1.0, 1.0, numFil, mxDOUBLE_CLASS,
                                    (int32_T)numFil, &d_emlrtRTEI, &st);
      for (loop_ub = 0; loop_ub < i2; loop_ub++) {
        b_st.site = &ob_emlrtRSI;
        sortNodes(&b_st, (real_T)loop_ub + 1.0, network->nodes,
                  crossings->filCross, asters->centers, astralNum, thisOrder);
        thisOrder_data = thisOrder->data;
        vectorUB = (int32_T)muDoubleScalarFloor(
            (((real_T)loop_ub + 1.0) - 1.0) / astralNum);
        if ((thisOrder->size[0] == 0) || (thisOrder->size[1] == 0)) {
          /*  0 nodes = filament doesn't touch any other asters! */
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &j_emlrtBCI,
                                          &st);
          }
          network->ends->data[loop_ub + network->ends->size[0]] = l;
          /*  treat unused filament as "right" dangling end */
        } else if (thisOrder->size[0] == 1) {
          real_T springLength;
          /*  1 node + astralNum >=2 -> spring between astral center & node */
          springCount++;
          if (vectorUB + 1 > numAsters_tmp) {
            emlrtDynamicBoundsCheckR2012b(vectorUB + 1, 1, numAsters_tmp,
                                          &u_emlrtBCI, &st);
          }
          if (((int32_T)thisOrder_data[0] < 1) ||
              ((int32_T)thisOrder_data[0] > maxdimlen)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)thisOrder_data[0], 1,
                                          maxdimlen, &t_emlrtBCI, &st);
          }
          numFil = network->nodes->data[(int32_T)thisOrder_data[0] - 1] -
                   asters->centers->data[vectorUB];
          b_numAsters[0] = numFil * numFil;
          numFil = network->nodes->data[((int32_T)thisOrder_data[0] +
                                         network->nodes->size[0]) -
                                        1] -
                   asters->centers->data[vectorUB + asters->centers->size[0]];
          b_numAsters[1] = numFil * numFil;
          springLength = sumColumnB(b_numAsters);
          b_st.site = &nb_emlrtRSI;
          if (springLength < 0.0) {
            emlrtErrorWithMessageIdR2018a(
                &b_st, &l_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
          }
          springLength = muDoubleScalarSqrt(springLength);
          thisOrder_data[0] += (real_T)c_numAsters;
          newSpring_size[0] = 1;
          newSpring_size[1] = 4;
          numFil = thisOrder_data[0];
          if (((int32_T)springCount < 1) ||
              ((int32_T)springCount > springs->size[0])) {
            emlrtDynamicBoundsCheckR2012b((int32_T)springCount, 1,
                                          springs->size[0], &m_emlrtBCI, &st);
          }
          a[0] = 1;
          a[1] = 4;
          emlrtSubAssignSizeCheckR2012b(&a[0], 2, &newSpring_size[0], 2,
                                        &emlrtECI, &st);
          springs_data[(int32_T)springCount - 1] = vectorUB + 1;
          springs_data[((int32_T)springCount + springs->size[0]) - 1] = numFil;
          springs_data[((int32_T)springCount + springs->size[0] * 2) - 1] =
              loop_ub + 1;
          springs_data[((int32_T)springCount + springs->size[0] * 3) - 1] =
              springLength;
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &emlrtBCI,
                                          &st);
          }
          network->ends->data[loop_ub + network->ends->size[0]] =
              l - springLength;
        } else {
          real_T springLength;
          real_T totSpringLength;
          /*  >=2 nodes = there are certainly springs along this filament */
          if (vectorUB + 1 > numAsters_tmp) {
            emlrtDynamicBoundsCheckR2012b(vectorUB + 1, 1, numAsters_tmp,
                                          &s_emlrtBCI, &st);
          }
          totSpringLength = 0.0;
          /*  like "Classical", but proximal segment to center is spring */
          /*  also need to account for augmented node indices */
          n = thisOrder->size[0] - 2;
          for (jdx = 0; jdx <= n; jdx++) {
            springCount++;
            if (jdx + 1 > thisOrder->size[0]) {
              emlrtDynamicBoundsCheckR2012b(jdx + 1, 1, thisOrder->size[0],
                                            &r_emlrtBCI, &st);
            }
            i = (int32_T)thisOrder_data[jdx];
            if ((i < 1) || (i > maxdimlen)) {
              emlrtDynamicBoundsCheckR2012b(i, 1, maxdimlen, &q_emlrtBCI, &st);
            }
            if (jdx + 2 > thisOrder->size[0]) {
              emlrtDynamicBoundsCheckR2012b(jdx + 2, 1, thisOrder->size[0],
                                            &p_emlrtBCI, &st);
            }
            i3 = (int32_T)thisOrder_data[jdx + 1];
            if ((i3 < 1) || (i3 > maxdimlen)) {
              emlrtDynamicBoundsCheckR2012b(i3, 1, maxdimlen, &o_emlrtBCI, &st);
            }
            numFil = network->nodes->data[i - 1] - network->nodes->data[i3 - 1];
            b_numAsters[0] = numFil * numFil;
            numFil = network->nodes->data[(i + network->nodes->size[0]) - 1] -
                     network->nodes->data[(i3 + network->nodes->size[0]) - 1];
            b_numAsters[1] = numFil * numFil;
            springLength = sumColumnB(b_numAsters);
            b_st.site = &lb_emlrtRSI;
            if (springLength < 0.0) {
              emlrtErrorWithMessageIdR2018a(
                  &b_st, &l_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                  "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
            }
            springLength = muDoubleScalarSqrt(springLength);
            numFil = thisOrder_data[jdx + 1];
            if ((!(thisOrder_data[jdx] > numFil)) &&
                ((!muDoubleScalarIsNaN(thisOrder_data[jdx])) ||
                 muDoubleScalarIsNaN(numFil))) {
              numFil = thisOrder_data[jdx];
            }
            numSpringsGuess = thisOrder_data[jdx + 1];
            if ((!(thisOrder_data[jdx] < numSpringsGuess)) &&
                ((!muDoubleScalarIsNaN(thisOrder_data[jdx])) ||
                 muDoubleScalarIsNaN(numSpringsGuess))) {
              numSpringsGuess = thisOrder_data[jdx];
            }
            if (((int32_T)springCount < 1) ||
                ((int32_T)springCount > springs->size[0])) {
              emlrtDynamicBoundsCheckR2012b(
                  (int32_T)springCount, 1, springs->size[0], &gb_emlrtBCI, &st);
            }
            springs_data[(int32_T)springCount - 1] =
                (real_T)c_numAsters + numFil;
            springs_data[((int32_T)springCount + springs->size[0]) - 1] =
                (real_T)c_numAsters + numSpringsGuess;
            springs_data[((int32_T)springCount + springs->size[0] * 2) - 1] =
                (real_T)loop_ub + 1.0;
            springs_data[((int32_T)springCount + springs->size[0] * 3) - 1] =
                springLength;
            totSpringLength += springLength;
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b(&st);
            }
          }
          /*  recording proximal segment as spring */
          springCount++;
          if (((int32_T)thisOrder_data[0] < 1) ||
              ((int32_T)thisOrder_data[0] > maxdimlen)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)thisOrder_data[0], 1,
                                          maxdimlen, &n_emlrtBCI, &st);
          }
          numFil = network->nodes->data[(int32_T)thisOrder_data[0] - 1] -
                   asters->centers->data[vectorUB];
          b_numAsters[0] = numFil * numFil;
          numFil = network->nodes->data[((int32_T)thisOrder_data[0] +
                                         network->nodes->size[0]) -
                                        1] -
                   asters->centers->data[vectorUB + asters->centers->size[0]];
          b_numAsters[1] = numFil * numFil;
          springLength = sumColumnB(b_numAsters);
          b_st.site = &kb_emlrtRSI;
          if (springLength < 0.0) {
            emlrtErrorWithMessageIdR2018a(
                &b_st, &l_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
          }
          springLength = muDoubleScalarSqrt(springLength);
          if (((int32_T)springCount < 1) ||
              ((int32_T)springCount > springs->size[0])) {
            emlrtDynamicBoundsCheckR2012b((int32_T)springCount, 1,
                                          springs->size[0], &eb_emlrtBCI, &st);
          }
          springs_data[(int32_T)springCount - 1] = (real_T)vectorUB + 1.0;
          springs_data[((int32_T)springCount + springs->size[0]) - 1] =
              (real_T)c_numAsters + thisOrder_data[0];
          springs_data[((int32_T)springCount + springs->size[0] * 2) - 1] =
              (real_T)loop_ub + 1.0;
          springs_data[((int32_T)springCount + springs->size[0] * 3) - 1] =
              springLength;
          totSpringLength += springLength;
          if ((int32_T)((uint32_T)loop_ub + 1U) > network->ends->size[0]) {
            emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 1U), 1,
                                          network->ends->size[0], &b_emlrtBCI,
                                          &st);
          }
          network->ends->data[loop_ub + network->ends->size[0]] =
              l - totSpringLength;
        }
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(&st);
        }
      }
    }
    emxFree_real_T(&st, &thisOrder);
    if ((int32_T)springCount == 0) {
      network->springs->size[0] = 0;
      network->springs->size[1] = 4;
    } else {
      if (springs->size[0] < 1) {
        emlrtDynamicBoundsCheckR2012b(1, 1, springs->size[0], &l_emlrtBCI, &st);
      }
      loop_ub = (int32_T)springCount;
      if (((int32_T)springCount < 1) ||
          ((int32_T)springCount > springs->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)springCount, 1, springs->size[0],
                                      &k_emlrtBCI, &st);
      }
      i = network->springs->size[0] * network->springs->size[1];
      network->springs->size[0] = (int32_T)springCount;
      network->springs->size[1] = 4;
      emxEnsureCapacity_real_T(&st, network->springs, i, &w_emlrtRTEI);
      for (i = 0; i < 4; i++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          network->springs->data[i2 + network->springs->size[0] * i] =
              springs_data[i2 + springs->size[0] * i];
        }
      }
    }
    emxFree_real_T(&st, &springs);
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (generateAstralNetwork.c) */
