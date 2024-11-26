#include <fstream>
#include <vector>
#include <iomanip>
#include "../lib/EcmData.h"
#include "../lib/generic_funcs.h"

typedef std::vector<double> dvec;
typedef std::vector<int> ivec;
typedef basics::vec2D<int> ivec2D;

namespace {
    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <input file>" << std::endl;
        return -1;
    }
    ivec2D qMomList {{0,0,0},{1,0,0},{1,1,0},{1,1,1},{2,0,0}};
}
int main(int argc, char** argv)
{
    if (argc != 2) {
        return usage(argc, argv);
    }
    std::string infile = argv[1];
    try {
        std::cout << "Reading in data" << std::endl;
        ecm::EcmData Data(infile);
        std::cout << "Finished reading in data\n" << std::endl;

        // Create the output file name as the input file name with "_out.dat" appended
        std::string outfile = infile;
        outfile.insert(outfile.find(".dat"), "_out_v2");
        // remove "inputs/" from the outfile name and replace with "data/"
        outfile.replace(0, 7, "data/");

        std::cout << "Now computing Q^2 values and writing to file: " << outfile << std::endl;
        Data.outputEvsQsq(outfile, qMomList);        
    }
    catch (std::string s) {
        std::cout << "\n\nERROR IN PROGRAM:\n" << s << std::endl;
    }
    return 0;
}