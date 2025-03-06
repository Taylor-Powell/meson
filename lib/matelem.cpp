#include <numbers>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "matelem.h"
#include "Wigner.h"


namespace matelem {
    state::state(std::string p, std::vector<int> threemom, double anis, int J, int P, int row, int Jz, bool current) {
        std::stringstream ss(p);
        ss >> V >> momstr >> irrep >> E >> Eerr;
        twopi_chiL = 2.0 * std::numbers::pi / (anis * V);
        mom = threemom;
        std::vector<double> threeMom = {mom[0] * twopi_chiL, mom[1] * twopi_chiL, mom[2] * twopi_chiL};
        fourMom = {E, threeMom[0], threeMom[1], threeMom[2]};
        double mom3_sq = basics::dot(threeMom, threeMom);
        params = p;
        spin = J;
        spinZ = Jz;
        parity = P;
        if (current) mState = 0.0;
        else mState = std::sqrt((basics::fourDot(fourMom, fourMom)).real()); 
        coeff = 1.0;
        etaTilde = parity * std::pow(-1, spin);
        irrepRow = row;
        sym = rotations::getSym(momstr);
        polVec = rotations::getPol4_Jz(E, mom3_sq, mom, Jz, sym, current);
    }

    void matelem::expandAllHelOps() {
        ExpandHelOps(init);
        ExpandHelOps(cur);
        ExpandHelOps(fin);
    }

    void matelem::subductAll(bool isHelState) {
        if (isHelState) {
            subductHelicityState(init);
            subductHelicityState(cur);
            subductHelicityState(fin);
        }
        else {
            subductJzState(init);
            subductJzState(cur);
            subductJzState(fin);
        }
    }
    
    void matelem::calcKinFactors() {
        std::vector<cd> kin;
        for (int i = 0; i < init.size(); i++) {
            for (int j = 0; j < cur.size(); j++) {
                for (int k = 0; k < fin.size(); k++) {
                    kin = kinFactors(init[i], cur[j], fin[k]);
                    #if 1
                    std::cout << "Kinematic factors for state " << i << " " << j << " " << k << ": ";
                    std::cout << "E1 = " << kin[0];
                    std::cout << ", C1 = " << kin[1] << std::endl;
                    #endif
                    kFactors.push_back(kin);
                }
            }
        }
    }

    cd matelem::getQsq(state& in, state& out) {
        std::vector<cd> qMom;
        std::cout << "qMom = ";
        for (int i = 0; i < 4; i++) {
            qMom.push_back(out.fourMom[i] - in.fourMom[i]);
            std::cout << out.fourMom[i] << " - " << in.fourMom[i] << " = " << qMom[i] << ", ";
        }
        std::cout << "qMom size = " << qMom.size() << std::endl;
        return (-1.0) * basics::fourDot(qMom, qMom);
    }

    void matelem::subductHelicityState(std::vector<state>& s) {
        if (s.size() != 1) throw std::string("subductState called with vector of size != 1.\n");
        if (s[0].etaTilde == 0) s[0].etaTilde = s[0].parity * std::pow(-1, s[0].spin);
        if (s[0].helicity == 0) s[0].coeff *= basics::subductHelicity(s[0].etaTilde, s[0].irrep, s[0].momstr, s[0].helicity, s[0].irrepRow);
        else {
            s[0].coeff *= basics::subductHelicity(s[0].etaTilde, s[0].irrep, s[0].momstr, s[0].helicity, s[0].irrepRow);
            s.push_back(s[0]);
            s[1].helicity = -s[1].helicity;
            s[1].coeff *= basics::subductHelicity(s[1].etaTilde, s[1].irrep, s[1].momstr, s[1].helicity, s[1].irrepRow);
        }
    }

    void matelem::subductJzState(std::vector<state>& s) {
        if (s.size() != 1) throw std::string("subductState called with vector of size != 1.\n");
        
    }
    
    void matelem::ExpandHelOps(std::vector<state>& s) {
        std::vector<state> s2;
        for (int i = 0; i < s.size(); i++) {
            state sTemp;
            std::vector<double> angles;
            std::vector<state> newStates;
            for (int i = 0; i <= 2 * s[i].spin + 1; i++) {
                sTemp = s[i];
                angles = rotations::getRotAngles(sTemp.momstr, sTemp.mom);
                sTemp.spinZ = -s[i].spin + i;
                sTemp.coeff *= std::conj(WignerD::Wigner_D(s[i].spin, s[i].spinZ, s[i].helicity, angles[0], angles[1], angles[2]));
                newStates.push_back(sTemp);
            }
            s2.insert(s2.end(), newStates.begin(), newStates.end());
        }
        s = s2;
    }

    cd matelem::getOmegaVal(state& in, state& out) {
        cd val = basics::fourDot(in.fourMom, out.fourMom);
        val -= std::pow(in.mState, 2) * std::pow(out.mState, 2);        
        return val;
    }

    std::vector<cd> matelem::kinFactors(state& in, state& cur, state& out) {
        std::vector<cd> kin;
        cd Omega = getOmegaVal(in, out);
        cd tempVal;

        // Calculate the E1 coefficient
        std::vector<cd> Ecoeff;
        for (int i = 0; i < 4; i++) {
            tempVal = basics::fourDot(in.fourMom, out.fourMom) * in.fourMom[i];
            tempVal -= std::pow(in.mState, 2) * out.fourMom[i];
            tempVal *= basics::fourDot(in.polVec, out.fourMom) / Omega;
            tempVal = out.polVec[i] - tempVal; 
            Ecoeff.push_back(tempVal);
        }
        kin.push_back(basics::fourDot(Ecoeff, cur.polVec));


        // Calculate the C1 coefficient
        std::vector<cd> Ccoeff;
        for (int i = 0; i < 4; i++) {
            tempVal = basics::fourDot(in.fourMom, out.fourMom);
            tempVal *= in.fourMom[i] * out.fourMom[i];
            tempVal -= std::pow(in.mState, 2) * out.fourMom[i];
            tempVal -= std::pow(out.mState, 2) * in.fourMom[i];
            tempVal *= basics::fourDot(in.polVec, out.fourMom) * in.mState / Omega;
            tempVal /= std::sqrt(basics::fourDot(cur.fourMom, cur.fourMom));
        }
        kin.push_back(basics::fourDot(Ccoeff, cur.polVec));

        return kin;
    }
}