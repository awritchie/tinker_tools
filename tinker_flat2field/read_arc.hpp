#ifndef tinker2field_read_arc_hpp
#define tinker2field_read_arc_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <limits>
#include "H5Cpp.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "tinker2field.hpp"

#define uindKfac (299.8391483043805 / 2.5852)  // converts from debye to uni    ts of KbT/(e- Angstrom)
#define frcKfac 1 // I haven't taken the time to figure this one out yet

class ReadArc
{
private:
    std::string name;
    int i1, i2;
    bool checkExt(std::string expected);
    std::vector<std::vector<float> > x1, x2;
    std::vector<float> a1, a2, mid, drop;
    bool isH5;
public:
    // Constructor
    ReadArc(std::string &fileName, int &a1, int &a2, std::string expected);
    // Deconstructor
    ~ReadArc();
    void ReadH5(std::string expected);
    int nFrames();
    std::vector<float> X(int i, int frame);
    float A1(int frame);
    float A2(int frame);
    float Mid(int frame);
    float Drop(int frame);
    void project(std::vector<std::vector<float> > &bondVector, int &frames, float k1, float k2, int NTHREADS);

    void fileHeader(FILE* file, std::string prefix);
    void writeLine(FILE* file, int linen);
    size_t skip2nextAtom(size_t current, std::vector<size_t> &ga1, std::vector<size_t> &ga2, bool smallStep);
     
};
    
#endif
