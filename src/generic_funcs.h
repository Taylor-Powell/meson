#ifndef __generic_funcs_h__
#define __generic_funcs_h__

#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <tuple>
#include <string>

namespace basics {
    
    /*
    * Set of momentum lists with allowed cubic rotations and their dims
    * Includes {000,001,011,111,002,012,112}
    */
    const int momList_001[][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    const int dimMomList_001 = 6;
    const int momList_002[][3] = {{2, 0, 0}, {-2, 0, 0}, {0, 2, 0}, {0, -2, 0}, {0, 0, 2}, {0, 0, -2}};
    const int dimMomList_002 = 6;
    const int momList_011[][3] = {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}, {1, -1, 0}, {0, 1, -1}, {-1, 0, 1}, {-1, 1, 0}, {0, -1, 1}, {1, 0, -1}, {-1, -1, 0}, {0, -1, -1}, {-1, 0, -1}};
    const int dimMomList_011 = 12;
    const int momList_111[][3] = {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, -1}};
    const int dimMomList_111 = 8;
    const int momList_012[][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}, {0, 1, -2}, {1, -2, 0}, {-2, 0, 1}, {0, -2, 1}, {-2, 1, 0}, {1, 0, -2}, {0, -1, 2}, {-1, 2, 0}, {2, 0, -1}, {0, 2, -1}, {2, -1, 0}, {-1, 0, 2}, {0, -1, -2}, {-1, -2, 0}, {-2, 0, -1}, {0, -2, -1}, {-2, -1, 0}, {-1, 0, -2}};
    const int dimMomList_012 = 24;
    const int momList_112[][3] = {{1, 1, 2}, {1, 2, 1}, {2, 1, 1}, {-1, 1, 2}, {-1, 2, 1}, {1, -1, 2}, {1, 2, -1}, {2, -1, 1}, {2, 1, -1}, {1, 1, -2}, {1, -2, 1}, {-2, 1, 1}, {-1, -1, 2}, {-1, 2, -1}, {2, -1, -1}, {-1, 1, -2}, {-1, -2, 1}, {1, -1, -2}, {1, -2, -1}, {-2, -1, 1}, {-2, 1, -1}, {-1, -1, -2}, {-1, -2, -1}, {-2, -1, -1}};
    const int dimMomList_112 = 24;

    /*
    * Given a momentum string in ascending order, returns the set of 
    * allowed permutations as a 2D vector by refernce
    */
    void momPerms(std::string mom, std::vector<std::vector<int>>& momList) {
        if (mom == "000") momList.push_back(std::vector<int>{0,0,0});
        else if (mom == "001") {
            for (int i = 0; i < dimMomList_001; i++) {
                momList.push_back(std::vector<int>{momList_001[i][0],momList_001[i][1],momList_001[i][2]});
            }
        }
        else if (mom == "011") {
            for (int i = 0; i < dimMomList_011; i++) {
                momList.push_back(std::vector<int>{momList_011[i][0],momList_011[i][1],momList_011[i][2]});
            }
        }
        else if (mom == "111") {
            for (int i = 0; i < dimMomList_111; i++) {
                momList.push_back(std::vector<int>{momList_111[i][0],momList_111[i][1],momList_111[i][2]});
            }
        }
        else if (mom == "002") {
            for (int i = 0; i < dimMomList_002; i++) {
                momList.push_back(std::vector<int>{momList_002[i][0],momList_002[i][1],momList_002[i][2]});
            }
        }
        else if (mom == "012") {
            for (int i = 0; i < dimMomList_012; i++) {
                momList.push_back(std::vector<int>{momList_012[i][0],momList_012[i][1],momList_012[i][2]});
            }
        }
        else if (mom == "112") {
            for (int i = 0; i < dimMomList_112; i++) {
                momList.push_back(std::vector<int>{momList_112[i][0],momList_112[i][1],momList_112[i][2]});
            }
        }
        // If not in allowed set, throw error.
        else {
            std::string errormsg = "Momentum " + mom + " not in allowed set {000,001,011,111,002,012,112}.\n";
            throw errormsg;
        }
        
    }
    
}


#endif