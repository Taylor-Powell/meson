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
#include <Eigen/Dense>

namespace {
    typedef std::complex<double> cd;
}

namespace basics {    
    template <typename T>
    using vec2D = std::vector<std::vector<T>>;

    // define a few macros
    #define MAX(x,y) (x>y ? x : y)
    #define MIN(x,y) (x<y ? x : y)

    template <typename T>
    inline std::string formatValue(T val) {
        static_assert(std::is_arithmetic<T>::value, "formatValue requires a numeric type.");
        std::string str = std::to_string(val);
        if (str.find(".") != std::string::npos) {
            str = str.substr(0, str.find(".") + 5);
            if ((str == "0.0000") || (str == "-0.0000")) return std::string("0");
        }
        return str;
    }

    template <typename T> // From Numerical Recipes
    T factorial(const int n) {
        static_assert(std::is_arithmetic<T>::value, "factorial requires a numeric type.");
        static int ntop = 4;
        static T a[33] = {1.0, 1.0, 2.0, 6.0, 24.0};
        int j;
        if (n < 0) throw std::domain_error("Negative factorial in routine factorial");
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
            throw std::invalid_argument("Attempting to dot vectors of unequal size.\n");
        for (int i = 0; i < v1.size(); i++) val += v1[i] * v2[i];
        return val;
    }

    template <typename T>
    inline T kDelta(const T i, const T j) {
        return (i == j) ? 1 : 0;
    }

    template <typename T>
    T fourDot(const std::vector<T>& v1, const std::vector<T>& v2) {
        if (v1.size() != 4 || v2.size() != 4)
            throw std::invalid_argument("Attempting fourDot vector(s) with size != 4.\n");
        return v1[0] * v2[0] - (v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]);
    }

    template <typename T>
    vec2D<T> matMult(const vec2D<T>& A, const vec2D<T>& B) {
        if (A.empty() || B.empty()) {
            throw std::invalid_argument("Matrices must be non-empty in basics::matMult.");
        }
        if (A[0].size() != B.size())
            throw std::invalid_argument("Attempting to multiply matrices of incompatible sizes.\n");
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
        if (A.empty()) {
            throw std::invalid_argument("Matrix must be non-empty in basics::matVecMult.");
        }
        if (A[0].size() != v.size())
            throw std::invalid_argument("Attempting to multiply matrix and vector of incompatible sizes.\n");
        // Initialize u as empty vector of size A.rows
        std::vector<T> u(A.size(), T());
        for (int i = 0; i < A.size(); i++)
            for (int j = 0; j < v.size(); j++)
                u[i] += A[i][j] * v[j];
        return u;
    }

    template <typename T>
    std::vector<T> rotVec(const std::vector<T>& vec, double phi, double theta, double psi) {
        if (vec.size() != 3) {
            throw std::invalid_argument("Vector size must be 3 in basics::rotVec.");
        }
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

    /** function to solve for the ordinary least-squares solution to an overdetermined system 
     * @param matrix: matrix of size (n, m) with n > m
     * @param b: vector of size n
     * @return: vector of size m with the least-squares solution
    */
    template <typename MatrixType, typename VectorType, typename T = typename MatrixType::Scalar>
    VectorType leastSquares(const MatrixType& matrix, const VectorType& b) {
        // If MatrixType is a vector of vectors, convert to Eigen::Matrix with scalar type <T>
        if constexpr (std::is_same_v<MatrixType, std::vector<std::vector<T>>>) {
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix(matrix.size(), matrix[0].size());
            for (size_t i = 0; i < matrix.size(); ++i) {
                for (size_t j = 0; j < matrix[i].size(); ++j) {
                    eigenMatrix(i, j) = matrix[i][j];
                }
            }
            return leastSquares(eigenMatrix, b); // Recursive call
        } 
        // If VectorType is a vector, convert to Eigen::Matrix with scalar type <T>
        else if constexpr (std::is_same_v<VectorType, std::vector<T>>) {
            Eigen::Matrix<T, Eigen::Dynamic, 1> eigenVector(b.size());
            for (size_t i = 0; i < b.size(); ++i) {
                eigenVector(i) = b[i];
            }
            Eigen::Matrix<T, Eigen::Dynamic, 1> result = leastSquares(matrix, eigenVector); // Recursive call

            // Convert the result back to std::vector<T>
            VectorType stdResult(result.size());
            for (size_t i = 0; i < result.size(); ++i) {
                stdResult[i] = result(i);
            }
            return stdResult;
        }

        // If matrix is square, solve directly
        if (matrix.rows() == matrix.cols()) {
            Eigen::Matrix<T, Eigen::Dynamic, 1> result = matrix.inverse() * b;

            // Convert the result back to VectorType if necessary
            if constexpr (std::is_same_v<VectorType, std::vector<T>>) {
                VectorType stdResult(result.size());
                for (size_t i = 0; i < result.size(); ++i) {
                    stdResult[i] = result(i);
                }
                return stdResult;
            }
            return result;
        }

        // Check if the matrix is underdetermined
        if (matrix.rows() < matrix.cols()) {
            throw std::invalid_argument("Matrix is underdetermined in basics::leastSquares.");
        }

        // Use QR decomposition for least-squares solution
        Eigen::Matrix<T, Eigen::Dynamic, 1> result = matrix.colPivHouseholderQr().solve(b);

        // Convert the result back to VectorType if necessary
        if constexpr (std::is_same_v<VectorType, std::vector<T>>) {
            VectorType stdResult(result.size());
            for (size_t i = 0; i < result.size(); ++i) {
                stdResult[i] = result(i);
            }
            return stdResult;
        }
        return result;
    }
    
    // Convenience function
    inline std::string getSym(const std::string momType) {
        if (momType == "00n") return "Dic4";
        else if (momType == "0nn") return "Dic2";
        else if (momType == "nnn") return "Dic3";
        else if (momType == "0mn") return "C40mn";
        else if (momType == "nnm") return "C4nnm";
        else if (momType == "000") return "OhD";
        else throw std::invalid_argument("Momentum " + momType + " not recognized in basics::getSym().\n");
    }

    /////////////////// Forward declarations ///////////////////

    // Overloaded function to get all permutations of a 3-momentum
    /** @param mom3_i: 3-momentum vector as 3-component vector, 0 <= x <= 9 */
    vec2D<int> getMomPerms(const std::vector<int>& mom3_i);
    /** @param mom3_i: 3-momentum vector as a string in form "xyz", 0 <= x <= 9 */
    vec2D<int> getMomPerms(const std::string mom3_i);

    // Overloaded function to get the irreps for a given momentum type
    std::vector<std::string> getIrreps(const std::vector<int> mom3_i, int parity, int spin);
    std::vector<std::string> getIrreps(const std::string momType, int parity, int spin);
    
    // Function to get the subduction coefficients for a given helicity
    double subductHelicity(int etaTilde, const std::string irrep, const std::string momType, int helicity, int irrepRow);

    // Simple function to check if momentum is <= 211 in all permutations */
    bool check3Mom(const std::vector<int> mom3_i);

    std::vector<int> getMom3_i(const std::string momStr);

    // Convenience function
    std::string getMomType(const std::vector<int> mom3_i);

    // Test function for Eigen
    void computeGramMatrix(const std::vector<std::vector<double>>& mat) {
        // Map mat to a MatrixXd
        Eigen::MatrixXd eigenMat(mat.size(), mat[0].size());
        for (size_t i = 0; i < mat.size(); ++i) {
            for (size_t j = 0; j < mat[i].size(); ++j) {
                eigenMat(i, j) = mat[i][j];
            }
        }
    }
}

#endif