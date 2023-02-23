#include "alloc.h"

void CheckLapackResult(int info, const char * routine) {
    if (info) {
        std::stringstream ss;
        ss << "Failure with LAPACK routine " << routine << ", info = " << info;
        throw std::runtime_error(ss.str());
    }
};

SexpWrapper::SexpWrapper(SEXP&& data) : data(PROTECT(data)) {

}

SexpWrapper::~SexpWrapper() {
    UNPROTECT_PTR(data);
}

SexpWrapper::operator SEXP&() {
    return data;
}
