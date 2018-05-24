#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void denCVBwEval(void *, void *, void *, void *, void *, void *);
extern void locCteSmoother(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void locCteWeights(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void locCuadSmoother(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void locLinSmoother(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void locLinWeights(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void locPolSmoother(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void locWeightsEvalxx(void *, void *, void *, void *, void *);
extern void looLocPolSmoother(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void parzenRossen(void *, void *, void *, void *, void *, void *, void *, void *);
extern void regCVBwEvalB(void *, void *, void *, void *, void *, void *, void *, void *);
extern void simpleSmoother(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void simpleSqSmoother(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"denCVBwEval",       (DL_FUNC) &denCVBwEval,        6},
    {"locCteSmoother",    (DL_FUNC) &locCteSmoother,    10},
    {"locCteWeights",     (DL_FUNC) &locCteWeights,      9},
    {"locCuadSmoother",   (DL_FUNC) &locCuadSmoother,   12},
    {"locLinSmoother",    (DL_FUNC) &locLinSmoother,    11},
    {"locLinWeights",     (DL_FUNC) &locLinWeights,      9},
    {"locPolSmoother",    (DL_FUNC) &locPolSmoother,    12},
    {"locWeightsEvalxx",  (DL_FUNC) &locWeightsEvalxx,   5},
    {"looLocPolSmoother", (DL_FUNC) &looLocPolSmoother, 10},
    {"parzenRossen",      (DL_FUNC) &parzenRossen,       8},
    {"regCVBwEvalB",      (DL_FUNC) &regCVBwEvalB,       8},
    {"simpleSmoother",    (DL_FUNC) &simpleSmoother,     9},
    {"simpleSqSmoother",  (DL_FUNC) &simpleSqSmoother,   8},
    {NULL, NULL, 0}
};

void R_init_locpol(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
