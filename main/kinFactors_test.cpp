#include <fstream>
#include <algorithm>
#include "../lib/outputKinFactors.h"

namespace {
    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <inputFile> <outputFile>" << std::endl;
        return -1;
    }
    typedef std::vector<double> dvec;
    typedef std::vector<int> ivec;
    typedef basics::vec2D<int> ivec2D;
    
    ivec2D qMomList {{0,0,0},{0,0,1},{0,1,1},{1,1,1},{0,0,2}};
}
int main(int argc, char** argv)
{
    #if 1
    // Testing readAndLoop class
    if (argc != 3) {
        return usage(argc, argv);
    }

    try {        
        kinFactors::Data d(argv[1]);
        // d.printParams();
        // d.outputKinematics(argv[2], qMomList);
    }
    catch (std::string s) {
        std::cout << "\n\nERROR IN PROGRAM:\n" << s << std::endl;
        return -1;
    }

    #endif

    return 0;
}