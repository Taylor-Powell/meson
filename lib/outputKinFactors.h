#ifndef __outputKinFactors_h__
#define __outputKinFactors_h__

#include <fstream>
#include <string>
#include <vector>
#include "generic_funcs.h"
#include "matelem.h"

namespace kinFactors {
    struct qTuple {
        std::vector<int> mom3_i;
        std::string irrep, momType;
        int irrepRow;
        double Q_sq;
    };
    
    class Data {
        public:
            // Constructors and Destructor
            Data() {}
            Data(std::string filename) { readData(filename); }
            ~Data() {}

            // Inline functions
            inline int getNumLevels() { return numLvls; };

            // Functions defined in read_Ecm_dat.cpp
            void readData(const std::string filename);
            void printParams();
            void outputKinematics(const basics::vec2D<int> qMomList);

        private:
            struct outState {
                int V;
                std::string momstr, momType, irrep, params, outfile;
                double E, Eerr;
                std::vector<int> mom3_i;
            };
            std::vector<outState> outStates;
            std::vector<std::string> outStrings;
            std::vector<std::vector<cd>> kFactors;
            double anis, at_mpi, at_inv;
            int numLvls, parity, spin, etaTilde;

            // Convenience functions for outputKinematics
            std::vector<qTuple> getqTuples(const basics::vec2D<int> qMomList);
            std::vector<std::pair<cd, int>> Data::getHelStates(int etaTilde, std::string momType, std::string irrep, int irrepRow, int spin);
            std::string getPiParamString(const std::vector<int> mom3_i, double anis, double at_mpi, double twopi_chiL, int V);
    };
}
#endif