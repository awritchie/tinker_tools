#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include "tinkerxyz.hpp"

#ifndef M_PI
#define M_PI (atan(1.d0)*4
#endif

typedef struct
{
    std::string outname;
    FILE *fp;

} analyze_t;



void analyze_frame_function(TinkerXYZ *frame, void *data) { 
    analyze_t *d = (analyze_t *)data;
    /* Open output file */
    std::ofstream xyz;
    xyz.open(d->outname+"_"+std::to_string(frame->nframes())+".arc");
    xyz << std::setw(6) << frame->natoms() << std::endl;
    if ((int)frame->box.size() == 6) {
        char buffer[1024];
        sprintf(buffer,"%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f",frame->box[0],frame->box[1],frame->box[2],frame->box[3],frame->box[4],frame->box[5]);
        xyz << buffer << std::endl;
    }
    for (int i=0; i<frame->natoms(); i++) {
        char buffer[1024];
        sprintf(buffer,"%6d  %-3s%12.6f%12.6f%12.6f%6d", i+1, frame->atom[i].c_str(), frame->xyz[i][0], frame->xyz[i][1], frame->xyz[i][2], frame->atype[i]);
        xyz << buffer;
        for (int j=0; j<(int)frame->bonded[i].size(); j++) {
            xyz << std::setw(6) << frame->bonded[i][j];
        }
        xyz << std::endl;
    }
    xyz.close();
}

int main(int argc, const char *argv[]) {
    /* Setup stuff, nothing to do with reading the tinker file */
    std::string usage = "\
    \n\
    Break a .arc file into individual frame files.\n\n\
    \t-x <coordinate file>\n\
    \t-o <out file name>\n\
    \t-b <first frame to read>\n\
    \t-e <last frame to read>\n\
    ";
    analyze_t d;
    
    /* Defaults */
    std::string filename;
    d.outname = "out";
    int begin = -1;
    int end = -1;
    if (argc>1) { 
        for (int i=1; i<argc; i++) { 
            if (strcmp("-h", argv[i]) == 0) { 
                std::cerr << usage << std::endl;
                exit(1);
            } else if (strcmp("-x", argv[i]) ==0) { 
                filename = argv[++i];
            } else if (strcmp("-o", argv[i]) ==0) { 
                d.outname = argv[++i];
            } else if (strcmp("-b", argv[i]) ==0) {
                begin = std::stoi(argv[++i]);
            } else if (strcmp("-e", argv[i]) ==0) { 
                end = std::stoi(argv[++i]);
            }
        }
    }
    /* Test beginning and end */
    if ((end < begin) && ((end != -1) && (begin != -1))) { 
        std::cerr << usage << std::endl;
        std::cerr << "\nERROR: The last frame (" << end << ") is before the first frame (" << begin << ").\n" << std::endl;
        exit(1);
    }
    if (begin < 0) { begin = 0; }

    /* Remove extension, if it exists, from outname */
    std::size_t m = d.outname.rfind(".");
    if (m != std::string::npos) { 
        d.outname = d.outname.substr(0,m);
    }

    /* Output the user selected options */
    std::cerr << "\n";
    std::cerr << "Reading " << filename << std::endl;
    std::cerr << "Reading from frame " << begin << " to frame ";
    if (end < 0) { std::cerr << "end\n";
    } else { std::cerr << end << std::endl; }
    std::cerr << "Writing to " << d.outname << ".xyz\n";
    std::cerr << std::endl;

    /* These few lines are the entirety of the tinker reading! */
    TinkerXYZ v(filename);
    v.isOctahedron();
    while (v.get_next_frame()) {
        if (v.nframes() >= begin) {
            analyze_frame_function(&v, &d);
        }
        if ((v.nframes() >= end) && (end > 0)) break;
    } v.close_tinker_read();
    std::cerr << std::endl;
    return 0;
}
