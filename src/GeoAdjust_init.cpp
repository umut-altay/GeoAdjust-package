#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void R_init_GeoAdjust(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

