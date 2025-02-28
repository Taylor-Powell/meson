#ifndef __readAndLoop_h__
#define __readAndLoop_h__

#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include "generic_funcs.h"
#include "matelem.h"

namespace readAndLoop {
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
            void readData(std::string filename);
            void printParams();
            void outputKinematics(std::string outfile, basics::vec2D<int> qMomList);

        private:
            std::vector<matelem::state> outStates;
            std::queue<matelem::matelem> matelems;
            double anis, at_mpi, at_mb1, at_inv;
            int numLvls, parity, spin, etaTilde;

            // Convenience functions for outputKinematics
            std::vector<qTuple> getqTuples(basics::vec2D<int> qMomList);
            std::string getPiParamString(std::vector<int> piMom, double anis, double at_mpi, double twopi_chiL);
    };
}
#endif