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

namespace basics {
    typedef std::complex<double> cd;
    cd j1 = cd(0,1);

    template <typename T>
    using vec2D = std::vector<std::vector<T>>;
    
    template <typename T>
    T dot(std::vector<T>& v1, std::vector<T>& v2) {
        T val = T();
        if (v1.size() != v2.size()) 
            throw std::string("Attempting to dot vectors of unequal size.\n");
        for (int i = 0; i < v1.size(); i++) val += v1[i] * v2[i];
        return val;
    }

    template <typename T>
    T kDelta(T i, T j) {
        return (i == j) ? 1 : 0;
    }

    template <typename T>
    T fourDot(std::vector<T>& v1, std::vector<T>& v2) {
        if (v1.size() != 4 || v2.size() != 4)
            throw std::string("Attempting fourDot vector(s) with size != 4.\n");
        return v1[0] * v2[0] - dot(v1, v2);
    }

    template <typename T>
    vec2D<T> matMult(vec2D<T>& A, vec2D<T>& B) {
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
    std::vector<T> matVecMult(vec2D<T>& A, std::vector<T>& v) {
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
    std::vector<T> rotVec (std::vector<T>& vec, double phi, double theta, double psi) {
        double cphi = std::cos(phi);
        double sphi = std::sin(phi);
        double ctheta = std::cos(theta);
        double stheta = std::sin(theta);
        double cpsi = std::cos(psi);
        double spsi = std::sin(psi);
        std::vector<std::vector<T>> rot1 = {{cphi,-sphi,0},{sphi,cphi,0},{0,0,1}};
        std::vector<std::vector<T>> rot2 = {{ctheta,0,-stheta},{0,1,0},{stheta,0,ctheta}};
        std::vector<std::vector<T>> rot3 = {{cpsi,-spsi,0},{spsi,cpsi,0},{0,0,1}};
        return matVecMult(matMult(matMult(rot1,rot2),rot3),vec);
    }

    /*
    * Function to search through vec2D list for a matching pair.
    *   - Returns false is a match is found
    *   - Default indices are {0,1}
    *   - Default precision is 1.0e-10
    */
    template <typename T>
    bool isUniquePair_vec2D(vec2D<T>& list, T v1, T v2, int i1=0, int i2=1, T eps = 1.0e-10) {
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
    int findMatchingPair_vec2D(vec2D<T>& list, T v1, T v2, int i1=0, int i2=1, T eps = 1.0e-10) {
        for (int i = 0; i < list.size(); i++)
            if ((std::abs(list[i][i1] - v1) < eps) && 
                (std::abs(list[i][i2] - v2) < eps))
                return i;
        return -1; 
    }

    // Forward declarations
    vec2D<int> getMomPerms(std::string mom);
    std::vector<std::string> getIrreps(int etaTilde, std::string mom, int helicity);
}


#endif