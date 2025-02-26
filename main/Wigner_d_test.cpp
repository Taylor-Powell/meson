#include <fstream>
#include "../lib/Wigner.h"

namespace {
    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <J> <M> <N> <beta>" << std::endl;
        return -1;
    }
}
int main(int argc, char** argv)
{
    #if 1
    // Testing Wigner_D function
    if (argc != 5) {
        return usage(argc, argv);
    }
    int J = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    int N = std::stoi(argv[3]);
    double beta = std::stod(argv[4]);
    std::cout << "Wigner_d(" << J << "," << M << "," << N << "," << beta << ") = " << WignerD::Wigner_d(J, M, N, beta) << std::endl;
    #endif

    return 0;
}