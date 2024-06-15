#include <R.h>
#include <Rinternals.h>

SEXP glue_series_lists(SEXP vexp) {
    size_t size, i, j, pnt;
    double* v;
    double* x;
    SEXP xexp;

    size = 0;

    // Rprintf("%d\n", length(vexp));

    for (i = 0; i < length(vexp); i++) {
        size += length(VECTOR_ELT(vexp, i));
    }

    // Rprintf("%d\n", size);

    xexp = PROTECT(allocVector(REALSXP, size));

    x = REAL(xexp);
    pnt = 0;
    for (i = 0; i < length(vexp); i++) {

        v = REAL(VECTOR_ELT(vexp, i));
        for (j = 0; j < length(VECTOR_ELT(vexp, i)); j++) {
            x[pnt++] = v[j];
        }
    }

    UNPROTECT(1);

    return xexp;
}
