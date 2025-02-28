#ifndef __wigner_h__
#define __wigner_h__

#include <complex>

namespace WignerD {
    double Wigner_d(int J, int M, int N, double beta);
    std::complex<double> Wigner_D(int J, int M, int N, double alpha, double beta, double gamma);
}

#endif