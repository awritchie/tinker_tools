#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "tinkerxyz.hpp"

#ifndef M_PI
#define M_PI (atan(1.d0)*4
#endif

typedef struct
{
    std::vector<std::vector<int> > dih;
    FILE *fp;
} analyze_t;


void analyze_frame_function(TinkerXYZ *frame, void *data) { 
    analyze_t *d = (analyze_t *)data;
    fprintf(d->fp,"%10i ",frame->nframes());
    for (int i=0; i<(int)d->dih.size(); i++) {
        fprintf(d->fp,"%8.2f ",frame->dihedral(d->dih[i][0],d->dih[i][1],d->dih[i][2],d->dih[i][3]));
    }
    fprintf(d->fp,"\n");
}

int main(int argc, const char *argv[]) {
    /* Setup stuff, nothing to do with reading the tinker file */
    std::string usage = "\
    This is just a template.\n\
    \n\
    It currently looked for hydrogen bonds and\n\
    calculated the radial distribution function\n\
    for the supplied trajectory file.\n\n\
    \t-x <coordinate file>\n\
    \t-n <index file containing dihedrals, separated with new lines>\n\
    \t-o <out file name>\n\
    \t-b <first frame to read>\n\
    \t-e <last frame to read>\n\
    ";
    analyze_t d;
    
    /* Defaults */
    std::string filename, ndxname;
    std::string outname = "out";
    int begin = -1;
    int end = -1;
    if (argc>1) { 
        for (int i=1; i<argc; i++) { 
            if (strcmp("-h", argv[i]) == 0) { 
                std::cerr << usage << std::endl;
                exit(1);
            } else if (strcmp("-x", argv[i]) ==0) { 
                filename = argv[++i];
            } else if (strcmp("-n", argv[i]) ==0) {
                ndxname = argv[++i];
            } else if (strcmp("-o", argv[i]) ==0) { 
                outname = argv[++i];
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
    std::size_t m = outname.rfind(".");
    if (m != std::string::npos) { 
        outname = outname.substr(0,m);
    }

    /* Output the user selected options */
    std::cerr << "\n";
    std::cerr << "Reading " << filename << std::endl;
    std::cerr << "Reading from frame " << begin << " to frame ";
    if (end < 0) { std::cerr << "end\n";
    } else { std::cerr << end << std::endl; }
    std::cerr << "Writing to " << outname << ".dih.\n";
    std::cerr << std::endl;
    
    /* Read index file */
    std::string line;
    std::ifstream ndxfile(ndxname);
    if (ndxfile.is_open()) {
        while ( getline(ndxfile,line)) {
            if (not line.empty()) {
                int index;
                std::vector<int> dih_ndx;
                std::stringstream linestream(line);
                while (linestream >> index) {
                    /* Decrement indices */
                    dih_ndx.push_back(index-1);
                }
                if ((int)dih_ndx.size() != 4) {
                    std::cerr << "\nError: Found " << dih_ndx.size() << " atoms for dihedral definition.\n";
                    std::cerr << line << std::endl;
                    exit(1);
                }
                d.dih.push_back(dih_ndx);
            }
        }
        ndxfile.close();
    } else {
        std::cerr << "Unable to open " << ndxname << std::endl;
        exit(1);
    }

    /* Open output file */
    d.fp = fopen((outname+".dih").c_str(), "w");

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
