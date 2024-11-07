// -*- C++ -*-
/* File: EvsQsq.h
 *
 * @brief This file contains 
 *
 */

#ifndef __EvsQsq_h__
#define __EvsQsq_h__

#include "read_Ecm_dat.h"
#include "generic_funcs.h"

namespace {
    typedef std::vector<int> ivec;
}

namespace ecm {
    double getQsq(ecm::energyLevel& l, ivec q, ivec mom);
    bool check3Mom(ivec mom3);
    void outputEvsQsqData(std::string outfile, EcmData& data, basics::vec2D<int> qMomList);
}

#endif