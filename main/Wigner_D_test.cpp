#include <fstream>
#include <iostream>
#include "../lib/Wigner.h"

namespace {
    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <J> <M> <N> <alpha> <beta> <gamma>" << std::endl;
        return -1;
    }
}
int main(int argc, char** argv)
{
    #if 1
    // Testing Wigner_D function
    if (argc != 7) {
        return usage(argc, argv);
    }
    int J = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    int N = std::stoi(argv[3]);
    double alpha = std::stod(argv[4]);
    double beta = std::stod(argv[5]);
    double gamma = std::stod(argv[6]);
    std::cout << "Wigner_D(" << J << "," << M << "," << N << "," << alpha << "," << beta << "," << gamma << ") = " << WignerD::Wigner_D(J, M, N, alpha, beta, gamma) << std::endl;
    #endif

    return 0;
}