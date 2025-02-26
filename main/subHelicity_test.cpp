#include <fstream>
#include "../lib/generic_funcs.h"

namespace {
    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <etaTilde> <irrep> <irrepRow> <mom> <helicity>" << std::endl;
        return -1;
    }
}
int main(int argc, char** argv)
{
    #if 1
    // Testing subductHelicity function
    if (argc != 6) {
        return usage(argc, argv);
    }
    int etaTilde = std::stoi(argv[1]);
    std::string irrep = argv[2];
    int irrepRow = std::stoi(argv[3]);
    std::string mom = argv[4];
    int helicity = std::stoi(argv[5]);
    std::cout << "SubductHelicity(" << etaTilde << ", " << irrep << ", " << mom << ", " << helicity << ", " << irrepRow << ") = " << basics::subductHelicity(etaTilde, irrep, mom, helicity, irrepRow) << std::endl;
    #endif

    return 0;
}