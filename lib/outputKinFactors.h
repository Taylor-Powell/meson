#ifndef __outputKinFactors_h__
#define __outputKinFactors_h__

#include <fstream>
#include <string>
#include "generic_funcs.h"
#include "matelem.h"

namespace kinFactors {
    struct qTuple {
        std::vector<int> qMom3;
        std::string irrep, momStr;
        int irrepRow, helicity;
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
            void outputKinematics(const std::string outfile, const basics::vec2D<int> qMomList);

        private:
            std::vector<matelem::state> outStates;
            double anis, at_mpi, at_inv;
            int numLvls, parity, spin, etaTilde;

            // Convenience functions for outputKinematics
            std::vector<qTuple> getqTuples(const basics::vec2D<int> qMomList);
            std::string getPiParamString(const std::vector<int> piMom, double anis, double at_mpi, double twopi_chiL, int V);
    };
}
#endif