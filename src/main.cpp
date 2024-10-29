#include <fstream>
#include "read_Ecm_dat.h"

namespace {
    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <dat filename>" << std::endl;
        return -1;
    }
}
int main(int argc, char** argv)
{    
    //////////////////////// Get E vs Qsq Plots ///////////////////////////
    #if 1
    if (argc != 2) {
        return usage(argc, argv);
    }
    std::string filename = argv[1];
    try {
        std::cout << "Reading in data\n";
        readEcm::EcmData E_structs(filename);
        std::cout << "Finished reading in data\n";


        readEcm::energyLevel l;
        while (!E_structs.isQueueEmpty()) {
            l = E_structs.popLevelInfo();
            std::cout << "mom = " << l.mom << std::endl;
        }
    }
    catch (std::string s) {
        std::cout << "\n\nERROR IN PROGRAM:\n" << s << std::endl;
    }

    #endif
    ///////////////////////////////////////////////////////////////////////
 

    return 0;
}