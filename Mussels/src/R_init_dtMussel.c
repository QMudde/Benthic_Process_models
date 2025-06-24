#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif


#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* register native routines ------------------------------------------------ */


/* NO .C calls */
/* NO .Call calls */

/* .Fortran calls */
void F77_NAME(initmuspar)  (void (* steadyparms)(int *, double *));
void F77_NAME(initmusforc) (void (* steadyforcs)(int *, double *));
void F77_NAME(musselmod)     (int *, double *, double *, double *, double *, int *);
void F77_NAME(musselcohort)  (int *, double *, double *, double *, double *, int *);

R_FortranMethodDef FEntries[] = {
  {"initmuspar",         (DL_FUNC) &F77_SUB(initmuspar),    1},
  {"initmusforc",        (DL_FUNC) &F77_SUB(initmusforc),  1},
  {"musselmod",          (DL_FUNC) &F77_SUB(musselmod),     6},
  {"musselcohort",       (DL_FUNC) &F77_SUB(musselcohort),     6},
  {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_dtMussel(DllInfo *dll) {

  R_registerRoutines(dll, NULL, NULL, FEntries, NULL);

  // the following line protects against accidentially finding entry points

  R_useDynamicSymbols(dll, FALSE); // disable dynamic searching
}
