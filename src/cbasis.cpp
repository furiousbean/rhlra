#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <fftw3.h>
#include "cbasis.h"
#include <complex>
#include "chorner.h"
#include "rotation_minimizer.h"
#include "calculate_basis.h"
#include "calculate_tangent_basis.h"
#include "calculate_pseudograd.h"
#include "calculate_sylvester_grad.h"

SEXP eval_basis(SEXP Nexp, SEXP GLRRexp, bool compensated) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;

    SEXP basis_exp = PROTECT(allocMatrix(CPLXSXP, N, r));
    SEXP unitroots_exp = PROTECT(allocVector(CPLXSXP, N));
    SEXP A_f_exp = PROTECT(allocVector(CPLXSXP, N));
    SEXP alpha_exp = PROTECT(allocVector(REALSXP, 1));
    SEXP answer_exp = PROTECT(allocVector(VECSXP, 5));

    SET_VECTOR_ELT(answer_exp, 0, basis_exp);
    SET_VECTOR_ELT(answer_exp, 1, unitroots_exp);
    SET_VECTOR_ELT(answer_exp, 2, A_f_exp);
    SET_VECTOR_ELT(answer_exp, 3, alpha_exp);
    SET_VECTOR_ELT(answer_exp, 4, GLRRexp);

    /* create names */
    SEXP names = PROTECT(allocVector(STRSXP, 5));
    SET_STRING_ELT(names, 0, mkChar("basis"));
    SET_STRING_ELT(names, 1, mkChar("unitroots"));
    SET_STRING_ELT(names, 2, mkChar("A_f"));
    SET_STRING_ELT(names, 3, mkChar("alpha"));
    SET_STRING_ELT(names, 4, mkChar("glrr"));

    /* assign names to list */
    setAttrib(answer_exp, R_NamesSymbol, names);

    std::complex<double>* basis_fourier = new std::complex<double>[N * r];

    if (compensated) {
        std::complex<double>* qrinvmat = new std::complex<double>[r * r];
        CalculateBasis<double, COMPENSATED_HORNER, ORTHOGONALIZATION> cb(N, r, REAL(GLRRexp),
            reinterpret_cast<std::complex<double>*>(COMPLEX(basis_exp)),
            basis_fourier,
            reinterpret_cast<std::complex<double>*>(COMPLEX(unitroots_exp)),
            reinterpret_cast<std::complex<double>*>(COMPLEX(A_f_exp)),
            REAL(alpha_exp),
            qrinvmat);
        cb.doWork();
        delete[] qrinvmat;
    } else {
        CalculateBasis<double> cb(N, r, REAL(GLRRexp),
            reinterpret_cast<std::complex<double>*>(COMPLEX(basis_exp)),
            basis_fourier,
            reinterpret_cast<std::complex<double>*>(COMPLEX(unitroots_exp)),
            reinterpret_cast<std::complex<double>*>(COMPLEX(A_f_exp)),
            REAL(alpha_exp),
            0);
        cb.doWork();
    }

    delete[] basis_fourier;

    UNPROTECT(6);

    return answer_exp;
}

SEXP eval_pseudograd(SEXP Nexp, SEXP GLRRexp,
    SEXP AFexp,
    SEXP ALPHAexp,
    SEXP TAUexp,
    SEXP SIGNALexp) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;
    int tau = INTEGER(TAUexp)[0] - 1;

    SEXP xexp = PROTECT(allocMatrix(CPLXSXP, N, r));

    CalculatePseudograd<double> cb(N, r, REAL(GLRRexp),
        reinterpret_cast<std::complex<double>*>(COMPLEX(AFexp)),
        reinterpret_cast<double*>(REAL(ALPHAexp)),
        tau,
        reinterpret_cast<double*>(REAL(SIGNALexp)),
        reinterpret_cast<std::complex<double>*>(COMPLEX(xexp)));

    cb.doWork();
    UNPROTECT(1);

    return xexp;
}

SEXP eval_sylvester_grad(SEXP Nexp, SEXP GLRRexp,
    SEXP AFexp,
    SEXP ALPHAexp,
    SEXP TAUexp,
    SEXP SIGNALexp) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;
    int tau = INTEGER(TAUexp)[0] - 1;

    SEXP xexp = PROTECT(allocVector(CPLXSXP, N));

    CalculateSylvesterGrad<double> cb(N, r, REAL(GLRRexp),
        reinterpret_cast<std::complex<double>*>(COMPLEX(AFexp)),
        reinterpret_cast<double*>(REAL(ALPHAexp)),
        tau,
        reinterpret_cast<double*>(REAL(SIGNALexp)),
        reinterpret_cast<std::complex<double>*>(COMPLEX(xexp)));

    cb.doWork();
    UNPROTECT(1);

    return xexp;
}

SEXP eval_tangent_basis(SEXP Nexp, SEXP GLRRexp, bool compensated) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;
    SEXP xexp = PROTECT(allocMatrix(CPLXSXP, N, 2 * r));

    if (compensated) {
        CalculateTangentBasis<double, COMPENSATED_HORNER> cb(N, r, REAL(GLRRexp),
            reinterpret_cast<std::complex<double>*>(COMPLEX(xexp)));
        cb.doWork();
    } else {
        CalculateTangentBasis<double> cb(N, r, REAL(GLRRexp),
            reinterpret_cast<std::complex<double>*>(COMPLEX(xexp)));
        cb.doWork();
    }

    UNPROTECT(1);

    return xexp;
}

extern "C" {
    SEXP eval_basis_compensatedC(SEXP Nexp, SEXP GLRRexp) {
        return eval_basis(Nexp, GLRRexp, true);
    }

    SEXP eval_basisC(SEXP Nexp, SEXP GLRRexp) {
        return eval_basis(Nexp, GLRRexp, false);
    }

    SEXP eval_pseudogradC(SEXP Nexp, SEXP GLRRexp,
        SEXP AFexp,
        SEXP ALPHAexp,
        SEXP TAUexp,
        SEXP SIGNALexp) {
        return eval_pseudograd(Nexp, GLRRexp,
        AFexp,
        ALPHAexp,
        TAUexp,
        SIGNALexp);
    }

    SEXP eval_sylvester_gradC(SEXP Nexp, SEXP GLRRexp,
        SEXP AFexp,
        SEXP ALPHAexp,
        SEXP TAUexp,
        SEXP SIGNALexp) {
        return eval_sylvester_grad(Nexp, GLRRexp,
        AFexp,
        ALPHAexp,
        TAUexp,
        SIGNALexp);
    }

    SEXP eval_tangent_basis_compensatedC(SEXP Nexp, SEXP GLRRexp) {
        return eval_tangent_basis(Nexp, GLRRexp, true);
    }

    SEXP eval_tangent_basisC(SEXP Nexp, SEXP GLRRexp) {
        return eval_tangent_basis(Nexp, GLRRexp, false);
    }
}
