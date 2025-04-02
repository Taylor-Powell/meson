#include <fstream>
#include <vector>
#include "generic_funcs.h"
#include "errorHandling.h"

namespace {
    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <v1> <v2> <v3> <v4>" << std::endl;
        return -1;
    }
}
int main(int argc, char** argv)
{
    #if 1
    // Testing subductHelicity function
    if (argc != 5) {
        return usage(argc, argv);
    }

    try {
        // Make a 2x2 matrix with v1 through v4 as the elements, which are double-valued
        std::vector<std::vector<double>> testMat = {
            {std::stod(argv[1]), std::stod(argv[2])},
            {std::stod(argv[3]), std::stod(argv[4])}
        };

        // Run it through basics::computeGramMatrix
        basics::computeGramMatrix(testMat);
    }   
    catch (...) {
        error::processException();
        return -1;
    }    
    #endif

    return 0;
}