#include <fstream>
#include <algorithm>
#include "outputKinFactors.h"
#include "errorHandling.h"
#include "timing.h"

namespace {
    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <inputFile>" << std::endl;
        return -1;
    }
    typedef std::vector<double> dvec;
    typedef std::vector<int> ivec;
    typedef basics::vec2D<int> ivec2D;
    
    ivec2D qMomList {{0,0,0},{0,0,1},{0,1,1},{1,1,1},{0,0,2}};
    // ivec2D qMomList {{0,0,1}};
}
int main(int argc, char** argv)
{
    if (argc != 2) {
        return usage(argc, argv);
    }

    try {
        timing::Timer timer;      
        kinFactors::Data d(argv[1]);
        d.outputKinematics(qMomList);
        timer.printElapsed("Execution time");
    }
    catch (...) {
        error::processException();
        return -1;
    }

    return 0;
}