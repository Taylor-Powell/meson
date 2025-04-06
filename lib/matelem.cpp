#include <numbers>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "matelem.h"
#include "Wigner.h"


namespace matelem {
    state::state(int V, std::string irrep, double E, double Eerr, std::vector<int> mom3_i, int J, int P, int row, int helicity, double twopi_chiL, bool current) 
        : V(V), irrep(irrep), E(E), Eerr(Eerr), mom3_i(mom3_i), spin(J), helicity(helicity), parity(P), irrepRow(row) {   
        if (mom3_i.size() != 3) {
            throw std::invalid_argument("Momentum vector must have exactly 3 components in state::state constructor.");
        }     
        momType = basics::getMomType(mom3_i);
        std::vector<double> mom3 = {mom3_i[0] * twopi_chiL, mom3_i[1] * twopi_chiL, mom3_i[2] * twopi_chiL};
        mom4 = {E, mom3[0], mom3[1], mom3[2]};
        double mom3_sq = basics::dot(mom3, mom3);
        if (current) mState = 0.0;
        else mState = std::sqrt((basics::fourDot(mom4, mom4)).real()); 
        coeff = 1.0;
        etaTilde = parity * std::pow(-1, spin);
        sym = basics::getSym(momType);
        polVec = rotations::getPol4_hel(E, mom3_sq, mom3_i, helicity, sym, current);
    }

    void matelem::subductAll(bool isHelState) {
        if (!isHelState) {
            throw std::invalid_argument("subductAll called with isHelState = false.\n");
        }

        // Helper lambda for debug output
        #if 0
        auto debugZeroCoeff = [&](const state& s, const std::string& type) {
            double eps = 1e-10;
            if (std::abs(s.coeff) < eps) {
                std::cout << "Zero for " << type << "(" << std::to_string(s.mom3_i[0]) 
                          << std::to_string(s.mom3_i[1]) << std::to_string(s.mom3_i[2]) << ") ";
                std::cout << s.irrep << "(" << s.irrepRow << ") etaTilde=" << s.etaTilde;
                std::cout << ", hel=" << s.helicity << std::endl;
            }
        };
        #else
        auto debugZeroCoeff = [&](const state& s, const std::string& type) {};
        #endif

        // Subduct helicity states and check for zero coefficients
        subductHelicityState(init);
        for (const auto& s : init) { debugZeroCoeff(s, "b1"); }

        subductHelicityState(cur);
        for (const auto& s : cur) { debugZeroCoeff(s, "current"); }

        subductHelicityState(fin);
        for (const auto& s : fin) { debugZeroCoeff(s, "pion"); }
    }
    
    bool matelem::calcKinFactors() {
        std::vector<cd> kin;
        std::vector<std::vector<int>> indices;
        std::vector<double> coeffs;
        double epsilon = 1e-10;
        int skipCount = 0, totalCount = 0;

        // Helper lambda to check if a subducted state has a zero coefficient
        auto zeroCoeff = [&](const state& s, const std::string& state) {
            if (std::abs(s.coeff) < epsilon) {
                #if 0
                std::cout << "Skipping " << state << " with irrep=" << s.irrep
                          << "(" << s.irrepRow << ") and momType=" << s.momType 
                          << std::endl;
                skipCount++;
                #endif
                return true;
            }
            return false;
        };

        for (int i = 0; i < init.size(); i++) {
            if (zeroCoeff(init[i], "b1")) continue;

            for (int j = 0; j < cur.size(); j++) {
                if (zeroCoeff(cur[j], "current")) continue;

                for (int k = 0; k < fin.size(); k++) {
                    if (zeroCoeff(fin[k], "pion")) continue;

                    coeffs.push_back(init[i].coeff * cur[j].coeff * fin[k].coeff);
                    kin.clear();
                    kin = kinFactors(init[i], cur[j], fin[k]);

                    if ((std::abs(kin[1]) < epsilon) && (std::abs(kin[2]) < epsilon)) {
                        #if 0
                        std::cout << "Skipping kinematic factors with zero coefficient." << std::endl;
                        skipCount++;
                        #endif

                        continue;
                    }

                    #if 0
                    std::cout << "Kinematic factors for state " << i << " " << j << " " << k << ": ";
                    std::cout << "E1 = " << kin[0];
                    std::cout << ", C1 = " << kin[1] << std::endl;
                    #endif

                    outStruct.kFactors.push_back(kin);
                    outStruct.helicity.push_back(cur[j].helicity);
                    indices.push_back({i, j, k});                    
                    totalCount++;
                }
            }
        }
        #if 0
        std::cout << "Skipped " << skipCount << " states with zero coefficient.    ";
        std::cout << "Calculated " << totalCount << " kinematic factors." << std::endl;
        #endif

        // If no kinematic factors were calculated, return false
        if (totalCount == 0) return false;

        // Format the output string
        // outStruct
        outStruct.outstring = formatOutString(coeffs, indices);
        return true;
    }

    cd matelem::getQsq(const state& in, const state& out) { // Checked
        if (in.mom4.size() != 4 || out.mom4.size() != 4) {
            throw std::invalid_argument("Four-momentum vectors must have exactly 4 components in matelem::getQsq().");
        }
        std::vector<cd> qMom;
        for (int i = 0; i < 4; i++) { qMom.push_back(out.mom4[i] - in.mom4[i]); }
        return (-1.0) * basics::fourDot(qMom, qMom);
    }

    void matelem::subductHelicityState(std::vector<state>& s) {
        if (s.size() != 1) throw std::invalid_argument("subductState called with vector of size != 1.\n");
        if (s[0].etaTilde == 0) s[0].etaTilde = s[0].parity * std::pow(-1, s[0].spin);
        if (s[0].helicity == 0) s[0].coeff *= basics::subductHelicity(s[0].etaTilde, s[0].irrep, s[0].momType, s[0].helicity, s[0].irrepRow);
        else {
            s[0].coeff *= basics::subductHelicity(s[0].etaTilde, s[0].irrep, s[0].momType, s[0].helicity, s[0].irrepRow);
            state sTemp = s[0];
            sTemp.helicity = -sTemp.helicity;
            sTemp.coeff *= basics::subductHelicity(sTemp.etaTilde, sTemp.irrep, sTemp.momType, sTemp.helicity, sTemp.irrepRow);
            s.push_back(sTemp);
        }
    }

    cd matelem::getOmegaVal(const state& in, const state& out) { // Checked
        if (in.mom4.size() != 4 || out.mom4.size() != 4) {
            throw std::invalid_argument("Four-momentum vectors must have exactly 4 components in matelem::getOmegaVal().");
        }
        cd val = std::pow(basics::fourDot(in.mom4, out.mom4), 2);
        val -= std::pow(in.mState, 2) * std::pow(out.mState, 2);
        // If Omega = 0, throw an error
        if (std::abs(val.real()) < 1e-10 && std::abs(val.imag()) < 1e-10) {
            // Throw message should include fourDot(in.mom4, out.mom4) and in.mState, out.mState
            throw std::runtime_error("Omega is zero in matelem::getOmegaVal(). FourDot(in.mom4, out.mom4) = " + std::to_string(basics::fourDot(in.mom4, out.mom4).real()) + " + " + std::to_string(basics::fourDot(in.mom4, out.mom4).imag()) + "i. m_in = " + std::to_string(in.mState) + ", m_out = " + std::to_string(out.mState));            
        }      
        return val;
    }

    std::vector<cd> matelem::kinFactors(const state& in, const state& cur, const state& out) { // Checked
        std::vector<cd> kin;
        cd Omega = getOmegaVal(in, out);
        cd tempVal;

        // Calculate Qsq
        kin.push_back(getQsq(in, out));

        // Lambda to check for NaN values
        auto checkNaN = [](const cd& val, const std::string& msg = "") {
            if (std::isnan(val.real()) || std::isnan(val.imag())) {
                throw std::runtime_error("Coefficient is NaN in matelem::kinFactors(). " + msg);
            }
        };

        // Calculate the E1 coefficient
        std::vector<cd> Ecoeff;
        for (int i = 0; i < 4; i++) { // Checked
            tempVal = basics::fourDot(in.mom4, out.mom4) * in.mom4[i];
            tempVal -= std::pow(in.mState, 2) * out.mom4[i];
            tempVal *= basics::fourDot(in.polVec, out.mom4) / Omega;
            tempVal = out.polVec[i] - tempVal;
            checkNaN(tempVal, "Computed Ecoeff[" + std::to_string(i) + "] for Ecoeff.");
            Ecoeff.push_back(tempVal);
        }
        kin.push_back(basics::fourDot(Ecoeff, cur.polVec));

        // Calculate the C1 coefficient
        std::vector<cd> Ccoeff;
        for (int i = 0; i < 4; i++) { // Checked
            tempVal = basics::fourDot(in.mom4, out.mom4);
            tempVal *= in.mom4[i] * out.mom4[i];
            tempVal -= std::pow(in.mState, 2) * out.mom4[i];
            tempVal -= std::pow(out.mState, 2) * in.mom4[i];
            tempVal *= basics::fourDot(in.polVec, out.mom4) * in.mState / Omega;
            tempVal /= std::sqrt(basics::fourDot(cur.mom4, cur.mom4));
            checkNaN(tempVal, "Computed Ccoeff[" + std::to_string(i) + "] for Ccoeff.");
            Ccoeff.push_back(tempVal);
        }
        kin.push_back(basics::fourDot(Ccoeff, cur.polVec));

        return kin;
    }

    // Helper function to format the output string
    std::string matelem::formatOutString(const std::vector<double>& coeffs, const std::vector<std::vector<int>>& indices) {
        std::ostringstream oss;

        // Format the header information
        oss << "Qsq=" << basics::formatValue(Qsq.real()) << ", ";
        oss << "b1(" << init[0].mom3_i[0] << init[0].mom3_i[1] << init[0].mom3_i[2] << ") ";
        oss << init[0].irrep << "(" << init[0].irrepRow << "), ";
        oss << "q(" << cur[0].mom3_i[0] << cur[0].mom3_i[1] << cur[0].mom3_i[2] << ") ";
        oss << cur[0].irrep << "(" << cur[0].irrepRow << "), ";
        oss << "pion(" << fin[0].mom3_i[0] << fin[0].mom3_i[1] << fin[0].mom3_i[2] << ") ";
        oss << fin[0].irrep << "(" << fin[0].irrepRow << ")\n";

        // Format the kinematic factors
        for (int i = 0; i < outStruct.kFactors.size(); i++) {
            oss << "     b1_hel=" << init[indices[i][0]].helicity << ", ";
            oss << "q_hel=" << cur[indices[i][1]].helicity << ", ";
            oss << "coeff=" << basics::formatValue(coeffs[i]) << ", ";
            oss << "E1=(" << basics::formatValue(outStruct.kFactors[i][1].real()) << ", " << basics::formatValue(outStruct.kFactors[i][0].imag()) << "), ";
            oss << "C1=(" << basics::formatValue(outStruct.kFactors[i][2].real()) << ", " << basics::formatValue(outStruct.kFactors[i][1].imag()) << ")";
            if (i != outStruct.kFactors.size() - 1) oss << "\n";
        }

        return oss.str();
    }
}

