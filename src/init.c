#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(log_sign)(double *lam2, double *lam3, double *gam, 
  int *nobs, int *nvars, double *x, double *y, int *jd, double *pf, 
  double *pf2, int *dfmax, int *pmax, int *nlam, double *flmin, double *ulam,
  double *eps, int *isd, int *maxit, int *nalam, double *b0, double *beta,
  int *ibeta, int *nbeta, double *alam, int *npass, int *jerr);

extern void F77_NAME(ls_sign)(double *lam2, double *lam3, double *gam, 
  int *nobs, int *nvars, double *x, double *y, int *jd, double *pf, 
  double *pf2, int *dfmax, int *pmax, int *nlam, double *flmin, double *ulam,
  double *eps, int *isd, int *maxit, int *nalam, double *b0, double *beta,
  int *ibeta, int *nbeta, double *alam, int *npass, int *jerr);

static const R_FortranMethodDef FortranEntries[] = {
  {"log_sign", (DL_FUNC) &F77_SUB(log_sign), 26},
  {"ls_sign", (DL_FUNC) &F77_SUB(ls_sign), 26},
  {NULL, NULL, 0}
};

void R_init_signpen(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
