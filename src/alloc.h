#ifndef ALLOC_H
#define ALLOC_H

#include <exception>
#include <sstream>

#include <fftw3.h>
#include <R.h>
#include <Rinternals.h>


void CheckLapackResult(int info, const char * routine);

template<class T> struct OrdinaryArrayDeleter {
    void operator()(T* p) {
        delete[] p;
    }
};

template<class T> struct FftwArrayDeleter {
    void operator()(T* p) {
        fftw_free(p);
    }
};

class SexpWrapper {
    SEXP data;

    public:
    SexpWrapper(SEXP&& data);

    ~SexpWrapper();

    operator SEXP&();
};

template<class T> T* FftwArrayAllocator(std::size_t N) {
    auto answer = reinterpret_cast<T*>(fftw_malloc(sizeof(T) * N));
    if (!answer) {
        std::stringstream ss;
        ss << "Failed to allocate memory for FFTW, size =  " << N;
        throw std::runtime_error(ss.str());
    }
    return answer;
};

#endif // ALLOC_H