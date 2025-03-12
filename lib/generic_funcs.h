#ifndef __generic_funcs_h__
#define __generic_funcs_h__

#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <tuple>
#include <string>
#include <cmath>
#include <complex>

namespace {
    typedef std::complex<double> cd;
}

namespace basics {    
    template <typename T>
    using vec2D = std::vector<std::vector<T>>;

    // define a few macros
    #define MAX(x,y) (x>y ? x : y)
    #define MIN(x,y) (x<y ? x : y)

    template <typename T> // From Numerical Recipes
    T factorial(const int n) {
        static int ntop = 4;
        static T a[33] = {1.0, 1.0, 2.0, 6.0, 24.0};
        int j;
        if (n < 0) throw std::string("Negative factorial in routine factorial");
        if (n > 32) return std::exp(std::lgamma(n + 1.0));
        while (ntop < n) {
            j = ntop++;
            a[ntop] = a[j] * ntop;
        }
        return a[n];
    }
    
    template <typename T>
    T dot(const std::vector<T>& v1, const std::vector<T>& v2) {
        T val = T();
        if (v1.size() != v2.size()) 
            throw std::string("Attempting to dot vectors of unequal size.\n");
        for (int i = 0; i < v1.size(); i++) val += v1[i] * v2[i];
        return val;
    }

    template <typename T>
    T kDelta(const T i, const T j) {
        return (i == j) ? 1 : 0;
    }

    template <typename T>
    T fourDot(const std::vector<T>& v1, const std::vector<T>& v2) {
        if (v1.size() != 4 || v2.size() != 4)
            throw std::string("Attempting fourDot vector(s) with size != 4.\n");
        return v1[0] * v2[0] - (v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]);
    }

    template <typename T>
    vec2D<T> matMult(const vec2D<T>& A, const vec2D<T>& B) {
        if (A[0].size() != B.size())
            throw std::string("Attempting to multiply matrices of incompatible sizes.\n");
        // Initialize C as empty matrix of size A.rows x B.cols
        vec2D<T> C(A.size(), std::vector<T>(B[0].size(), T()));
        for (int i = 0; i < A.size(); i++)
            for (int j = 0; j < B[0].size(); j++)
                for (int k = 0; k < B.size(); k++)
                    C[i][j] += A[i][k] * B[k][j];
        return C;
    }

    template <typename T>
    std::vector<T> matVecMult(const vec2D<T>& A, const std::vector<T>& v) {
        if (A[0].size() != v.size())
            throw std::string("Attempting to multiply matrix and vector of incompatible sizes.\n");
        // Initialize u as empty vector of size A.rows
        std::vector<T> u(A.size(), T());
        for (int i = 0; i < A.size(); i++)
            for (int j = 0; j < v.size(); j++)
                u[i] += A[i][j] * v[j];
        return u;
    }

    template <typename T>
    std::vector<T> rotVec (const std::vector<T>& vec, double phi, double theta, double psi) {
        double cphi = std::cos(phi);
        double sphi = std::sin(phi);
        double ctheta = std::cos(theta);
        double stheta = std::sin(theta);
        double cpsi = std::cos(psi);
        double spsi = std::sin(psi);
        vec2D<T> rot1 = {{cphi,-sphi,0},{sphi,cphi,0},{0,0,1}};
        vec2D<T> rot2 = {{ctheta,0,-stheta},{0,1,0},{stheta,0,ctheta}};
        vec2D<T> rot3 = {{cpsi,-spsi,0},{spsi,cpsi,0},{0,0,1}};
        vec2D<T> rot = matMult(rot1, rot2);
        rot = matMult(rot, rot3);
        return matVecMult(rot, vec);
    }

    /*
    * Function to search through vec2D list for a matching pair.
    *   - Returns false is a match is found
    *   - Default indices are {0,1}
    *   - Default precision is 1.0e-10
    */
    template <typename T>
    bool isUniquePair_vec2D(const vec2D<T>& list, T v1, T v2, int i1=0, int i2=1, T eps = 1.0e-10) {
        for (int i = 0; i < list.size(); i++)
            if ((std::abs(list[i][i1] - v1) < eps) && 
                (std::abs(list[i][i2] - v2) < eps))
                return false;
        return true;
    }
    /*
    * Function to search through vec2D list for a matching pair. 
    *   - v1 and v2 are the values to find matches for
    *   - i1 and i2 are the column indices to search over
    *   - Returns -1 if no match is found
    *   - otherwise returns row index of FIRST match
    *   - Default indices are {0,1}
    *   - Default precision is 1.0e-10
    */
    template <typename T>
    int findMatchingPair_vec2D(const vec2D<T>& list, T v1, T v2, int i1=0, int i2=1, T eps = 1.0e-10) {
        for (int i = 0; i < list.size(); i++)
            if ((std::abs(list[i][i1] - v1) < eps) && 
                (std::abs(list[i][i2] - v2) < eps))
                return i;
        return -1; 
    }

    // Forward declarations
    vec2D<int> getMomPerms(const std::string mom);
    std::vector<std::string> getIrreps(int etaTilde, const std::string mom, int helicity);
    double subductHelicity(int etaTilde, const std::string irrep, const std::string mom, int helicity, int irrepRow);
    /** Simple function to check if momentum is <= 211 in all permutations */
    bool check3Mom(const std::vector<int> mom3);
}


#endif