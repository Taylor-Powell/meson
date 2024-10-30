// -*- C++ -*-
/* File: EvsQsq.h
 *
 * Description: 
 *
 * This file contains 
 *
 */

#ifndef __EvsQsq_h__
#define __EvsQsq_h__

#include <cmath>
#include "read_Ecm_dat.h"
#include "generic_funcs.h"

namespace {
    typedef std::vector<int> ivec;
}

namespace ecm {
    double getQsq(ecm::energyLevel& l, ivec q, ivec mom) {
        double Eval = sqrt(pow(l.at_mpi,2) + pow(l.twopi_chiL, 2) * basics::dot(mom, mom));
        return -pow(l.E - Eval, 2) + pow(l.twopi_chiL, 2) * basics::dot(q, q);
    }

    bool check3Mom(ivec mom3) {
        if (((abs(mom3[0]) < 3) && (abs(mom3[0]) < 2) && (abs(mom3[0]) < 2))
            || ((abs(mom3[0]) < 2) && (abs(mom3[0]) < 3) && (abs(mom3[0]) < 2)) 
            || ((abs(mom3[0]) < 2) && (abs(mom3[0]) < 2) && (abs(mom3[0]) < 3)))
            return true;
        return false;
    }

    void outputEvsQsqData(std::string outfile, EcmData& data, basics::vec2D<int> qMomList) {
        std::ofstream out (outfile);
        double val;
        energyLevel l;
        basics::vec2D<int> pMomList;
        ivec qmom, pmom, mom3;
        std::tuple<double,int,double,int> pts(0.0, 0, 0.0, 1);
        basics::vec2D<double> outvals;
        ivec outcount;
        int index;

        if (!out) {
            std::string errormsg = "Failed to open file in ";
            errormsg += __func__;
            throw errormsg;
        }
        while (!data.isQueueEmpty()) {
            l = data.popLevelInfo();
            pMomList = basics::momPerms(l.mom);
            for (int i = 0; i < qMomList.size(); i++) {
                qmom = qMomList[0];
                for (int j = 0; j < pMomList.size(); j++) {
                    pmom = pMomList[0];
                    mom3 = { pmom[0]-qmom[0], pmom[1]-qmom[1], pmom[2]-qmom[2]};
                    if (!check3Mom(mom3)) continue;
                    val = getQsq(l, qmom, mom3);
                    index = basics::findMatchingPair_vec2D(outvals, l.E, val, 0,1);
                    if (index == -1) outvals.push_back({l.E, val});
                    else outcount[i] += 1;
                }
            }
        }
        

    }
}

#endif