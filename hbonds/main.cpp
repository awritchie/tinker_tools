#include <iostream>
#include <cmath>
#include <string>
#include "tinkerxyz.hpp"

#ifndef M_PI
#define M_PI (atan(1.d0)*4
#endif

typedef struct
{
    float sum, cutoff, width, maxdist;
    int ci, ni, hi, nbins;
    int owtype;
    std::vector<float> r, cnh, nho;
    std::vector<int> rdf_hist;
    FILE *fp, *rdf;

} analyze_t;

float dist(const std::vector<float> &a, const std::vector<float> &b) {
    float dx = a[0] - b[0];
    float dy = a[1] - b[1];
    float dz = a[2] - b[2];
    float rr = dx*dx + dy*dy + dz*dz;
    return std::sqrt(rr);
}

float dot(const std::vector<float> &a, const std::vector<float> &b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

float ang(const std::vector<float> &a, const std::vector<float> &b, const std::vector<float> &c) { 
    /*
                 c
                /
               /
        a-----b
    */
    std::vector<float> ba = { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
    std::vector<float> bc = { c[0] - b[0], c[1] - b[1], c[2] - b[2] };
    float rba = std::sqrt(dot(ba,ba));
    float rbc = std::sqrt(dot(bc,bc));
    for (int i=0; i<3; i++) {
        ba[i] /= rba;
        bc[i] /= rbc;
    }
    return std::acos(dot(ba,bc))*180/M_PI;
}

void analyze_frame_function(TinkerXYZ *frame, void *data) { 
    analyze_t *d = (analyze_t *)data;
    std::vector<float> cd = frame->xyz[d->ci];
    std::vector<float> ne = frame->xyz[d->ni];
    for (int i=0; i<frame->natoms(); i++) {
        if (frame->atype[i] == d->owtype) { 
            float oh1r = frame->distance(i+1, d->ni);
            float oh2r = frame->distance(i+2, d->ni);
            float owr = std::min(oh1r,oh2r);
            int bin0 = owr/d->width +1;
            if (bin0 < d->nbins) { 
                d->rdf_hist[bin0]++;
                d->sum++;
                }
            if (owr < d->cutoff+2) { 
                float r1 = frame->distance(i+1, d->ni);
                float r2 = frame->distance(i+2, d->ni);
                //r1 = dist(frame->xyz[i+1], ne);
                //r2 = dist(frame->xyz[i+1], ne);
                if ( (r1<d->cutoff) || (r2<d->cutoff) ) { 
                    if (r1<r2) { 
                        d->r.push_back(r1);
                        d->cnh.push_back(frame->angle(d->ci, d->ni, i+1));
                        d->nho.push_back(frame->angle(d->ni, i+1, i));
                    } else {
                        d->r.push_back(r2);
                        d->cnh.push_back(frame->angle(d->ci, d->ni, i+2));
                        d->nho.push_back(frame->angle(d->ni, i+2, i));
                    }
                    d->hi++;
                    if (d->fp) {
                        fprintf(d->fp, "%10i %10i %12.4f %12.4f %12.4f\n", frame->nframes(), i, d->r[d->hi-1], d->cnh[d->hi-1], d->nho[d->hi-1]);
                    }
                }
            }
            i+=2;
        }
    }
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
    \t-cutoff <maximum hbonding distance>\n\
    \t-ni <index of cnc nitrogen\n\
    \t-ci <index of cnc carbon\n\
    \t-ow <atom type of water oxygen\n\
    \t-width <bin width for RDF>\n\
    \t-max_dist <maximum distance for RDF>\n\
    \t-o <out file name>\n\
    \t-b <first frame to read>\n\
    \t-e <last frame to read>\n\
    ";
    analyze_t d;
    
    /* Defaults */
    std::string filename;
    std::string outname = "out";
    d.maxdist = 20.;
    d.width = 0.01;
    d.owtype = 247;
    d.cutoff = 2.45;
    int begin = -1;
    int end = -1;
    if (argc>1) { 
        for (int i=1; i<argc; i++) { 
            if (strcmp("-h", argv[i]) == 0) { 
                std::cerr << usage << std::endl;
                exit(1);
            } else if (strcmp("-x", argv[i]) ==0) { 
                filename = argv[++i];
            } else if (strcmp("-cutoff", argv[i])==0) { 
                d.cutoff = std::stof(argv[++i]);
            } else if (strcmp("-ni", argv[i]) ==0) {
                d.ni = std::stoi(argv[++i]);
            } else if (strcmp("-ci", argv[i]) ==0) {
                d.ci = std::stoi(argv[++i]);
            } else if (strcmp("-ow", argv[i]) ==0) { 
                d.owtype = std::stoi(argv[++i]);
            } else if (strcmp("-width", argv[i]) ==0) { 
                d.width = std::stof(argv[++i]);
            } else if (strcmp("-max_dist", argv[i]) ==0) {
                d.maxdist = std::stof(argv[++i]);
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
    std::cerr << "CNC CD, NE: " << d.ci << ", " << d.ni << std::endl;
    std::cerr << "Water Oxygen Type: " << d.owtype << std::endl;
    std::cerr << "Maximum NE-HW distance: " << d.cutoff << " Angstoms.\n";
    std::cerr << "Creating RDF up to " << d.maxdist << " Angstroms.\n";
    std::cerr << "Reading from frame " << begin << " to frame ";
    if (end < 0) { std::cerr << "end\n";
    } else { std::cerr << end << std::endl; }
    std::cerr << "Writing to " << outname << ".hb and " << outname << ".rdf.\n";
    std::cerr << std::endl;
    
    /* Initialize variables */
    d.sum = 0, d.hi =0;
    d.nbins = d.maxdist / d.width + 1;
    d.rdf_hist = std::vector<int> (d.nbins, 0);
    /* Decrement indices */
    d.ni--; d.ci--;
    if (d.ni < 0 || d.ci < 0) { 
        std::cerr << "\nERROR: -ci and -ni have not been assigned.\n" << std::endl;
        exit(1);
    } else if (d.ni == d.ci) { 
        std::cerr << "\nERROR: -ci (" << d.ci << ") and -ni (" << d.ni << ") should not have the same index.\n" << std::endl;
        exit(1);
    }

    /* Open output file */
    d.fp = fopen((outname+".hb").c_str(), "w");

    /* These few lines are the entirety of the tinker reading! */
    TinkerXYZ v(filename);
    v.isOctahedron();
    while (v.get_next_frame()) {
        if (v.nframes() >= begin) {
            analyze_frame_function(&v, &d);
        }
        if ((v.nframes() >= end) && (end > 0)) break;
    } v.close_tinker_read();

    /* Finishing writing output of analysis */
    std::cout << "\nThere are " << d.hi << " hbonds and " << v.nframes() << " frames.\n";
    fclose(d.fp);
    /* RDF stuff now */
    d.rdf = fopen((outname+".rdf").c_str(), "w");
    float factor = (4./3.) * M_PI * v.nframes() * d.sum / v.volume();
    for (int i=0; i<(int)d.rdf_hist.size(); i++) {
        float rupper = (i+1)*d.width;
        float rlower = rupper - d.width;
        float expect = factor * (pow(rupper,3) - pow(rlower,3));
        fprintf(d.rdf, "%.2f %18.6e\n", rupper, d.rdf_hist[i]/expect);
    }
    fclose(d.rdf);
    return 0;
}
