#ifndef tinker2field_read_prm_hpp
#define tinker2field_read_prm_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include "H5Cpp.h"

class ReadPrm
{
private:
    std::vector<float> alpha;
    std::vector<float> q;
    std::vector<int> atom;
    std::vector<int> type;
    void getTypes(std::string &crdName);
    void getPrm(std::string &prmName);
public:
    // Constructor
    ReadPrm(std::string &prmName, std::string &crdName, int &a1, int &a2);
    // Deconstructor
    ~ReadPrm();
    int Atom(int i);
    int Type(int i);
    float Alpha(int i);
    float Q(int i);
};



#endif
