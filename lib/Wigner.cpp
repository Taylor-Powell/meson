#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include "Wigner.h"
#include "generic_funcs.h"

namespace {
    auto fact = basics::factorial<double>;
}

namespace WignerD {
    /** Calculates the Wigner d-functions. 
     *
     * NOTE also: the args j1, m1, etc. are 2*spin,
     *
     * \returns \f$d^{j}_{m,n}(\beta)\f$.
     *
     */
    double Wigner_d(int J, int M, int N, double beta) {
        int temp_M, k, k_low, k_hi;
        double const_term = 0.0, sum_term = 0.0, d = 1.0;
        int m_p_n, j_p_m, j_p_n, j_m_m, j_m_n;
        int kmn1, kmn2, jmnk, jmk, jnk;
        double kk;

        if (J < 0 || std::abs (M) > J || std::abs (N) > J) {
            std::string errmsg = "\nWignerD::Wigner_d - You have entered an illegal number for J, M, N.";
            errmsg += "\nMust follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J.";
            errmsg += "\nJ = " + std::to_string(J) +  " M = " + std::to_string(M) + " N = " + std::to_string(N);
            throw std::string(errmsg);
            return 0.0;
        }
    
        if (beta < 0) {
            beta = std::fabs (beta);
            temp_M = M;
            M = N;
            N = temp_M;
        }

        m_p_n = (M + N) / 2;
        j_p_m = (J + M) / 2;
        j_m_m = (J - M) / 2;
        j_p_n = (J + N) / 2;
        j_m_n = (J - N) / 2;
    
        kk = fact(j_p_m) * fact(j_m_m) * fact(j_p_n) * fact(j_m_n);
        const_term = std::pow(-1.0, j_p_m) * std::sqrt(kk);  
    
        k_low = MAX(0, m_p_n);
        k_hi = MIN(j_p_m, j_p_n);

        for (k = k_low; k <= k_hi; k++) {
            kmn1 = 2 * k - (M + N) / 2;
            jmnk = J + (M + N) / 2 - 2 * k;
            jmk = (J + M) / 2 - k;
            jnk = (J + N) / 2 - k;
            kmn2 = k - (M + N) / 2;
                
            sum_term += std::pow(-1.0, k) *
            ((std::pow (std::cos(beta / 2.0), kmn1)) * (std::pow(std::sin(beta / 2.0), jmnk))) /
            (fact (k) * fact (jmk) * fact (jnk) * fact (kmn2));
        }

        d = const_term * sum_term;
        return d;
    }

    /** Calculates the Wigner D-functions. 
     *
     * NOTE also: the args j1, m1, etc. are 2*spin,
     *
     * \returns \f$D^{j}_{m,n}(\alpha, \beta, \gamma)\f$.
     *
     */
    std::complex<double> Wigner_D(int J, int M, int N, double alpha, double beta, double gamma) {
        const std::complex<double> j1(0.0, 1.0);

        if (J < 0 || abs (M) > J || abs (N) > J) {
            std::string errmsg = "\nWignerD::Wigner_D - you have entered an illegal number for J, M, N.";
            errmsg += "\nMust follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J.";
            errmsg += "\nJ = " + std::to_string(J) +  " M = " + std::to_string(M) + " N = " + std::to_string(N);
            throw std::string(errmsg);
            return 0.0;
        }

        double expPhase = - alpha * double(M) / 2.0 - gamma * double(N) / 2.0;
        auto D = std::cos(expPhase) + j1 * std::sin(expPhase);
        D = D * Wigner_d(J, M, N, beta);

        return D;
    }
}