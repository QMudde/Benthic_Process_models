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
void F77_NAME(initmacroalgae) (void (* funparms)(int *, double *));
void F77_NAME(initforcm)      (void (* funforcs)(int *, double *));
void F77_NAME(macroalgae)     (int *, double *, double *, double *, double *, int *);

R_FortranMethodDef FEntries[] = {
  {"initmacroalgae",     (DL_FUNC) &F77_SUB(initmacroalgae), 1},
  {"initforcm",          (DL_FUNC) &F77_SUB(initforcm),      1},
  {"macroalgae",         (DL_FUNC) &F77_SUB(macroalgae),     6},
  {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_dtPosidonia(DllInfo *dll) {

  R_registerRoutines(dll, NULL, NULL, FEntries, NULL);

  // the following line protects against accidentially finding entry points

  R_useDynamicSymbols(dll, FALSE); // disable dynamic searching
}
