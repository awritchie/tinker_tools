#ifndef tinkerxyz_hpp
#define tinkerxyz_hpp

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <cstring>
#include <vector>
#include <sstream>
#include <limits>
#include <cmath>

#ifndef M_PI
#define M_PI atan(1)/4
#endif

#define RADIAN (180./M_PI)

class TinkerXYZ
{
private:
    int numatoms, numframes;
    FILE *file;
    std::string filename;
    bool hasBox;
    bool isBox(char *c);
    std::vector<std::vector<float> > _xyz;
    std::vector<std::string> _atom;
    std::vector<int> _atype;
    std::vector<std::vector<int> > _bonded;
    std::vector<float> _box;
    std::vector<int> threshold, increment;

    bool octahedron, triclinic, monoclinic, orthogonal;
    float vdot(const std::vector<float> &a, const std::vector<float> &b);
    std::vector<float> vcross(const std::vector<float> &a, const std::vector<float> &b);
public:
    TinkerXYZ(const std::string &file);
    ~TinkerXYZ();
    void open_tinker_read();
    void close_tinker_read();
    void read_tinker_structure();
    void isOctahedron();
    std::vector<float> pair_minimum_vector(const int &a, const int &b);
    std::vector<float> pair_minimum_vector(const std::vector<float> &a, const std::vector<float> &b);
    float distance(const int &a, const int &b);
    float distance(const std::vector<float> &a, const std::vector<float> &b);
    float angle(const int &a, const int &b, const int &c);
    float angle(const std::vector<float> &a, const std::vector<float> &b, const std::vector<float> c);
    float dihedral(const int &a, const int &b, const int &c, const int &d);
    float dihedral(const std::vector<float> &a, const std::vector<float> &b, const std::vector<float> &c, const std::vector<float> &d);
    float volume();
    bool get_next_frame();
    int natoms();
    int nframes();

    const std::vector<std::vector<float> > &xyz;
    const std::vector<std::string> &atom;
    const std::vector<int> &atype;
    const std::vector<std::vector<int> > &bonded;
    const std::vector<float> &box;
};

#endif
