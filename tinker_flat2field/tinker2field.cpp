#include "tinker2field.hpp"

int main(int argc, const char * argv[]){
    std::string usage = "\
\tThis program reads in tinker coordinate files along with force files\n\
and/or induced dipole files and prints out the electrostatic field along a bond\n\
vector given by two atoms.\n\n\
\tThe field will be projected along the bond vector pointing from atom 1\n\
to atom 2.  When an induced dipole moment file is provided, the field is\n\
calculated as E_i = (induced dipole_i)/alpha_i, where alpha_i is the\n\
polarizability on atom i. When a force file is provided, the field is\n\
calculated as E_i = (force_i)/q_i, where q_i is the charge on atom i.\n\n\
Usage: tinker2field\n\
\t-x <coordinate file>\n\
\t-p <parameter file>\n\
\t-a1 <atom 1 index>\n\
\t-a2 <atom 2 index>\n\
\t-u <induced dipole file>\n\
\t-f <force file>\n\
\t-nt <number of threads>\n\
\t-o <out file: default out.txt>\n\
    ";
    
    // Read arguments
    std::string uindName, frcName, crdName, prmName, sa1, sa2;
    std::string snt = "0", out="out.txt";
    int NTHREADS, minFrames;
    bool doUind = false, doFrc = false;
    if (argc > 1) {
        for (int i=1; i<argc; i++) {
            if (strcmp("-h", argv[i]) == 0) {
                std::cerr << usage << std::endl;
                exit(1);
            }
            else if (strcmp("-x", argv[i]) ==0) {
                crdName = argv[++i];
            }
            else if (strcmp("-p", argv[i]) ==0) {
                prmName = argv[++i];
            }
            else if (strcmp("-u", argv[i]) ==0) {
                //if (argv[i+1])
                uindName = argv[++i];
                doUind = true;
            }
            else if (strcmp("-f", argv[i]) ==0) {
                frcName = argv[++i];
                doFrc = true;
            }
            else if (strcmp("-a1", argv[i]) ==0) {
                sa1 = argv[++i];
            }
            else if (strcmp("-a2", argv[i]) ==0) {
                sa2 = argv[++i];
            }
            else if (strcmp("-nt", argv[i]) ==0) {
                snt = argv[++i];
            }
            else if (strcmp("-o", argv[i]) ==0) {
                out = argv[++i];
            }
        }
    }

    int minNarg = 10, maxNarg = 16;
    if ((argc < minNarg+1) || (argc > maxNarg+1)) {
        std::cerr << usage << std::endl;
        std::cerr << "\nError: Unexpected number of arguments.\n\t" << argc -1 << " arguments recieved, expected " << minNarg << " or " << maxNarg << ".\n" << std::endl;
        exit(1);
    }
    if ( not doUind && not doFrc ) {
        std::cerr << "\nError: No induced dipole moment or force files specified.\n";
        exit(1);
    }
    
    // Cast the atoms into ints
    int a1 = stoi(sa1);
    int a2 = stoi(sa2);
    if (a1 == 0) {
        std::cerr << "\nError: Could not convert -a1 to integer > 0.\n\n";
        exit(1);
    }
    if (a2 == 0) {
        std::cerr << "\nError: Could not convert -a2 to integer > 0.\n\n";
        exit(1);
    }
    if (a1 == a2) {
        std::cerr << "\nError: -a1 and -a2 are the same.\n\n";
        exit(1);
    }
    
#ifdef _OPENMP
    if (snt == "0") {
        NTHREADS = omp_get_max_threads();
    }
    else {
        NTHREADS = stoi(snt);
    }
    std::cout << "\nUsing OMP with up to " << NTHREADS << " threads.\n";
    omp_set_num_threads(NTHREADS);
#endif
    
    // Print out information about what is being examined.
    std::cout << "\nReading coordinates from < " << crdName << " >\n";
    std::cout << "Reading parameters from < " << prmName << " >\n";
    std::cout << "Looking at the bond vector pointing from atoms < " << a1 << " to " << a2 << " >\n";
    if (doUind) {
        std::cout << "Reading induced dipoles from < " << uindName << " >\n"; }
    if (doFrc) {
        std::cout << "Reading forces from < " << frcName << " >\n"; }

    ReadPrm* prmDat;
    ReadArc* crdDat;
    ReadArc* uindDat;
    ReadArc* frcDat;

#ifdef _OPENMP
    int sTHREADS = 2;
    if (doUind) { sTHREADS++; }
    if (doFrc) {sTHREADS++;}
    sTHREADS = std::min(sTHREADS,NTHREADS);
    std::cout << "\nUsing " << sTHREADS << " threads for file reading.\n";
//#endif
#pragma omp parallel sections num_threads(sTHREADS)
{
#pragma omp section
#endif
    {
        // Read parameter data
        prmDat = new ReadPrm(prmName, crdName, a1, a2);
    }
#ifdef _OPENMP
#pragma omp section
#endif
    {
        // Read coordinate data
        crdDat = new ReadArc(crdName, a1, a2, "coordinate");
    }
#ifdef _OPENMP
#pragma omp section
#endif
    {
        // Read induced dipole data
        if (doUind) {
            uindDat = new ReadArc(uindName, a1, a2, "induced dipole");
        }
    }
#ifdef _OPENMP
#pragma omp section
#endif
    {
        // Read force data
        if (doFrc) {
            frcDat = new ReadArc(frcName, a1, a2, "force");
        }
    }
#ifdef _OPENMP
}
#endif

    minFrames = crdDat->nFrames();
    // Make vector of normalized bond vectors
    std::vector<std::vector<float> > bondVector (crdDat->nFrames(),std::vector<float> (3,0));
#ifdef _OPENMP
#pragma omp parallel for num_threads(NTHREADS)
#endif
    for (int i=0; i<crdDat->nFrames(); i++) {
        float r2 = 0;
        for (int j=0; j<3; j++) {
            bondVector[i][j] = crdDat->X(2,i)[j] - crdDat->X(1,i)[j];
            r2 += bondVector[i][j] * bondVector[i][j];
        }
        float rr = 1./std::sqrt(r2);
        for (int j=0; j<3; j++) {
            bondVector[i][j] *= rr;
        }
    }
    
    // Project the induced dipole moments onto the normalized bond vectors and
    // divide by alpha
    if (doUind) {
        uindDat->project(bondVector, minFrames, uindKfac/prmDat->Alpha(0), uindKfac/prmDat->Alpha(1), NTHREADS);
    }
    if (doFrc) {
        frcDat->project(bondVector, minFrames, frcKfac/prmDat->Q(0), frcKfac/prmDat->Q(1), NTHREADS);
    }
    // Write to file 
    if (doUind or doFrc) {
        FILE* file = fopen(out.c_str(), "w");
        std::cout << "Writing output to " << out << ".\n";
        fprintf(file, ";%9s","frame");
        if (doUind) { uindDat->fileHeader(file,"u"); }
        if (doFrc) { frcDat->fileHeader(file,"f"); }
        fprintf(file,"\n");
        for (int i=0; i<minFrames; i++){
            fprintf(file,"%10i",i+1);
            if (doUind) { uindDat->writeLine(file,i); }
            if (doFrc) { frcDat->writeLine(file,i); }
            fprintf(file,"\n");
        }
        fclose(file);
    }
    
    // Clean up
    delete prmDat;
    delete crdDat;
    delete uindDat;
    delete frcDat;
    return 0;
}
