#include <fstream>
#include <vector>
#include <iomanip>
#include "../lib/read_Ecm_dat.h"
#include "../lib/generic_funcs.h"
#include "../lib/EvsQsq.h"

typedef std::vector<double> dvec;
typedef std::vector<int> ivec;
typedef basics::vec2D<int> ivec2D;

namespace {
    int usage1(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <input file>" << std::endl;
        return -1;
    }
    int usage2(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <input file> <output file>" << std::endl;
        return -1;
    }
    ivec2D qMomList {{0,0,0},{1,0,0},{1,1,0},{1,1,1},{2,0,0}};
}
int main(int argc, char** argv)
{    
    ///////////////// Test reading in data from file /////////////////////
    #if 1
    if (argc != 2) {
        return usage1(argc, argv);
    }
    std::string infile = argv[1];
    try {
        std::cout << "Reading in data\n";
        ecm::EcmData E_structs(infile);
        std::cout << "Finished reading in data\n\n";
        
        E_structs.printParams();
    }
    catch (std::string s) {
        std::cout << "\n\nERROR IN PROGRAM:\n" << s << std::endl;
    }
    #endif
    ///////////////////////////////////////////////////////////////////////


    //////////////////// Test function momPerms() /////////////////////////
    #if 0
    if (argc != 2) {
        return usage1(argc, argv);
    }
    std::string infile = argv[1];
    try {
        std::cout << "Reading in data\n";
        ecm::EcmData E_structs(infile);
        std::cout << "Finished reading in data\n\n";
       
        std::string str("011");
        basics::vec2D<int> momList = basics::momPerms(str);

        std::cout << "momList for (mom = \"" << str << "\") is:" << std::endl;
        for (int i = 0; i < momList.size(); i++) {
            for (int j = 0; j < momList[i].size(); j++) {
                std::cout << momList[i][j] << "  ";
            }
            std::cout << std::endl;
        }
    }
    catch (std::string s) {
        std::cout << "\n\nERROR IN PROGRAM:\n" << s << std::endl;
    }
    #endif
    ///////////////////////////////////////////////////////////////////////


    ////////////////// Test matchingPair() functions //////////////////////
    #if 0
    if (argc != 2) {
        return usage1(argc, argv);
    }
    std::string infile = argv[1];
    try {
        std::cout << "Reading in data\n";
        ecm::EcmData E_structs(infile);
        std::cout << "Finished reading in data\n\n";
        
        ecm::energyLevel l;
        int row = 0;
        bool b;
        basics::vec2D<double> data(E_structs.getNumLevels(), dvec(2));
        while (!E_structs.isQueueEmpty()) {
            l = E_structs.popLevelInfo();
            data[row][0] = l.E;
            data[row][1] = l.err;   
            row++;         
        }
        // Test findMatchingPair_vec2D()
        row = basics::findMatchingPair_vec2D(data, data[0][0], data[0][1]);
        if (row == -1) std::cout << "No match found.\n";
        else std::cout << "Match found in row " << row << std::endl;
        row = basics::findMatchingPair_vec2D(data, data[32][0], data[32][1]);
        if (row == -1) std::cout << "No match found.\n";
        else std::cout << "Match found in row " << row << std::endl;
        row = basics::findMatchingPair_vec2D(data, data[15][0], data[16][1]);
        if (row == -1) std::cout << "No match found.\n";
        else std::cout << "Match found in row " << row << std::endl;

        // Test isUniquePair_vec2D()
        b = basics::isUniquePair_vec2D(data, data[0][0], data[0][1]);
        if (b) std::cout << "Pair is unique." << std::endl;
        if (!b) std::cout << "Pair is not unique." << std::endl;
        b = basics::isUniquePair_vec2D(data, data[32][0], data[32][1]);
        if (b) std::cout << "Pair is unique." << std::endl;
        if (!b) std::cout << "Pair is not unique." << std::endl;
        b = basics::isUniquePair_vec2D(data, data[15][0], data[16][1]);
        if (b) std::cout << "Pair is unique." << std::endl;
        if (!b) std::cout << "Pair is not unique." << std::endl;


    }
    catch (std::string s) {
        std::cout << "\n\nERROR IN PROGRAM:\n" << s << std::endl;
    }
    #endif
    ///////////////////////////////////////////////////////////////////////


    ////////////////// Test matchingPair() functions //////////////////////
    #if 0
    if (argc != 3) {
        return usage2(argc, argv);
    }
    std::string infile = argv[1];
    std::string outfile = argv[2];
    try {
        std::cout << "Reading in data\n";
        ecm::EcmData E_structs(infile);
        std::cout << "Finished reading in data\n\n";
    }
    catch (std::string s) {
        std::cout << "\n\nERROR IN PROGRAM:\n" << s << std::endl;
    }
    #endif
    ///////////////////////////////////////////////////////////////////////
 

    return 0;
}