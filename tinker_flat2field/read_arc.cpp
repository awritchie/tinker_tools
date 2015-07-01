#include "read_arc.hpp"

ReadArc::ReadArc(std::string &fileName, int &a1, int &a2, std::string expected)
{
    name = fileName;
    i1 = a1;
    i2 = a2;
    if (not checkExt(expected)) {
        exit(1);
    }
    if (isH5) {
        ReadH5(expected);
    }
    else {
       std::vector<size_t> pos_a1 (2,0);
       std::vector<size_t> pos_a2 (2,0);
       int linen = 0, natoms = 0;
       int skip = abs(a1-a2);
       std::string line;
       std::ifstream file(fileName);
       if (file.is_open()) {
           // The fist line has the number of atoms 
           getline(file,line);
           natoms = stoi(line);
           size_t position = file.tellg();
    
           linen++;
           while (std::getline(file,line)) { 
               if (not line.empty())
               {
                   std::stringstream linestream(line);
                   std::vector<std::string> lineVector;
                   std::string column;
                   while ( linestream >> column ) {
                       lineVector.push_back(column);
                   }
                   if (lineVector.size() > 2) {
                       std::vector<float> xyz (3,0);
                       if (stod(lineVector[0]) == stoi(lineVector[0])) {
                           // Check atom 1
                           if (stoi(lineVector[0]) == a1) {
                               for (int i=0; i<3; i++) {
                                   xyz[i] = stod(lineVector[i+2]);
                               }
                               if (x1.size() == 0) pos_a1[0] = position;
                               else if (x1.size() == 1) pos_a1[1] = position;
                               x1.push_back(xyz);
                               size_t skipTo = skip2nextAtom(position, pos_a1, pos_a2, (a1 < a2));
                               if (skipTo > 0) file.seekg(skipTo, file.beg);
                           }
                           // Check atom 2
                           if (stoi(lineVector[0]) == a2) {
                               for (int i=0; i<3; i++) {
                                   xyz[i] = stod(lineVector[i+2]);
                               }
                               if (x2.size() == 0) pos_a2[0] = position;
                               else if (x2.size() == 1) pos_a2[1] = position;
                               x2.push_back(xyz);
                               size_t skipTo = skip2nextAtom(position, pos_a1, pos_a2, (a2 < a1));
                               if (skipTo > 0) file.seekg(skipTo, file.beg);
                           }
                       }
                   }
               }
               linen++;
               position = file.tellg();
           }
       }
       else {
           fprintf(stderr,"\nError: Cannot open %s.\n",name.c_str());
           exit(1);
       }
    }
    nFrames();
    fprintf(stderr, "%s has %d frames\n", name.c_str(), nFrames());
}

ReadArc::~ReadArc(void){}

void ReadArc::ReadH5(std::string expected) {
    const H5std_string FILE_NAME(name);
    H5::H5File *file;
    try
    {
        H5::Exception::dontPrint();
        file = new H5::H5File(FILE_NAME, H5F_ACC_RDWR);
    }
    // catch failure caused by the H5File operations
    catch( H5::FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( H5::DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( H5::DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( H5::DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    
    std::string X;
    if (strcmp(expected.c_str(),"coordinate") == 0 ) { X = "X"; }
    else if (strcmp(expected.c_str(),"force") == 0 ) { X = "F"; }
    else if (strcmp(expected.c_str(),"induced dipole") == 0 ) { X = "U"; }
    else if (strcmp(expected.c_str(),"velocity") == 0 ) { X = "V"; }

    H5::Group *grp = new H5::Group(file->openGroup(X));
    int attr_nframes[1];
    H5::Attribute attr = grp->openAttribute("nframes");
    attr.read(H5::PredType::NATIVE_INT, attr_nframes);
    int nframes = attr_nframes[0];
    
    x1 = std::vector<std::vector<float> > (nframes, std::vector<float> (3,0));
    x2 = std::vector<std::vector<float> > (nframes, std::vector<float> (3,0));
    
    H5::DataSet *dataset;
    H5::DataSpace *dataspace;
    hsize_t a1_offset[2] = { (hsize_t)i1-1, 0 }; // offset by -1 b/c .xyz numbering starts
    hsize_t a2_offset[2] = { (hsize_t)i2-1, 0 }; // at 1 while array indexing starts at 0
    hsize_t count[2] = { 1, 3 };
    std::vector<float> atom1 (3,0), atom2 (3,0);
    const int one = 1;
    hsize_t row_dims[1] = { 3 };
    H5::DataSpace mspace2(one, row_dims);
    for (int i=0; i<nframes; i++) {
        try {
            dataset = new H5::DataSet(grp->openDataSet(std::to_string(i+1)));
            dataspace = new H5::DataSpace(dataset->getSpace());
            dataspace->selectHyperslab(H5S_SELECT_SET, count, a1_offset);
            dataset->read(&atom1[0], H5::PredType::NATIVE_FLOAT, mspace2, *dataspace);
            dataspace->selectHyperslab(H5S_SELECT_SET, count, a2_offset);
            dataset->read(&atom2[0], H5::PredType::NATIVE_FLOAT, mspace2, *dataspace);
            x1[i] = atom1;
            x2[i] = atom2;
            dataset->close();
            delete dataset;
            delete dataspace;
        }
        // catch failure caused by the H5File operations
        catch( H5::FileIException error )
        {
            error.printError();
            std::exit(1);
        }
        // catch failure caused by the DataSet operations
        catch( H5::DataSetIException error )
        {
            error.printError();
            std::exit(1);
        }
        // catch failure caused by the DataSpace operations
        catch( H5::DataSpaceIException error )
        {
            error.printError();
            std::exit(1);
        }
        // catch failure caused by the DataSpace operations
        catch( H5::DataTypeIException error )
        {
            error.printError();
            std::exit(1);
        }
        // catch failure caused by the Group operations
        catch( H5::GroupIException error )
        {
            error.printError();
            std::exit(1);
        }
    }
    file->close();
    return;
}
int ReadArc::nFrames() {
    if (x1.size() == x2.size()) {
        return (int)x1.size();
    }
    else {
        fprintf(stderr,"\nError: The number of frames for atom 1 do not match the number of frames for atom 2 in %s.\n",name.c_str());
        exit(1);
    }
}

std::vector<float> ReadArc::X(int i, int frame) {
    if (i==1) { return x1[frame]; }
    else if (i==2) { return x2[frame]; }
    else { 
        fprintf(stderr,"Error: Expected to look at atom 1 or atom 2, not atom %d.\n");
        exit(1);
    }
}

float ReadArc::A1(int frame) {
    if (frame < a1.size()) { return a1[frame]; }
    else { 
        fprintf(stderr,"\nError: Frame requested (%d) exceeds number of existing frames (%d).\n",frame,a1.size());
        exit(1);
    }
}

float ReadArc::A2(int frame) {
    if (frame < a2.size()) { return a2[frame]; }
    else { 
        fprintf(stderr,"\nError: Frame requested (%d) exceeds number of existing frames (%d).\n",frame,a2.size());
        exit(1);
    }
}
float ReadArc::Mid(int frame) {
    if (frame < mid.size()) { return mid[frame]; }
    else { 
        fprintf(stderr,"\nError: Frame requested (%d) exceeds number of existing frames (%d).\n",frame,mid.size());
        exit(1);
    }
}
float ReadArc::Drop(int frame) {
    if (frame < drop.size()) { return drop[frame]; }
    else { 
        fprintf(stderr,"\nError: Frame requested (%d) exceeds number of existing frames (%d).\n",frame,drop.size());
        exit(1);
    }
}
bool ReadArc::checkExt(std::string expected) {
    isH5 = false;
    int extStart = name.rfind(".") + 1;
    int ext_ = name.rfind("_");
    int length = name.size();
    std::string ext = name.substr(extStart,name.size() - extStart);
    std::string ext_u = ext;
    if (ext_ > 0) {
        ext_u = name.substr(extStart,ext_ - extStart);
    };
    const char* lastChar = &ext_u.c_str()[ext_u.size()-1];
    if (strcmp("coordinate",expected.c_str()) == 0) {
        if (strcmp("xyz", ext_u.c_str()) == 0) { return true; }
        if (strcmp("arc", ext_u.c_str()) == 0) { return true; }
        bool num = true;
        for (int i=0; i<ext.size(); i++) {
            if (std::isalpha(ext.c_str()[i])) {
                num=false;
                break;
            }
        }
        if (num) { return num; }
    }
    else if (strcmp("induced dipole",expected.c_str()) == 0) {
        if (strcmp("uind", ext_u.c_str()) == 0) { return true; }
        if (strcmp("u", lastChar) == 0) { return true; }
    }
    else if (strcmp("force",expected.c_str()) == 0) {
        if (strncmp("frc", ext.c_str(),3) == 0) { return true; }
        if (strcmp("f", lastChar) == 0) { return true; }
    }
    if (strcmp("h5", ext.c_str()) == 0) {
        isH5 = true;
        return true;
    }
    fprintf(stderr,"\nError: Unexpected extention (.%s) for a %s file\n.",ext.c_str(),expected.c_str());
    return false;
}

void ReadArc::project(std::vector<std::vector<float> > &bondVector, int &frames, float k1, float k2, int NTHREADS = 1) {
    int nframes = x1.size();
    if (frames != nframes) {
        fprintf(stderr, "\n Warning!\tThe number of frames in %s (%d) does not\n\t\tmatch the number of frames in the coordinate dataset (%d).\n\t\tAssuming the smaller of the two is the most complete.\n", name.c_str(), nframes, (int)bondVector.size());
        frames = std::min(frames,nframes);
    }
    a1   = std::vector<float> (frames,0);
    a2   = std::vector<float> (frames,0);
    mid = std::vector<float> (frames,0);
    drop = std::vector<float> (frames,0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(NTHREADS)
#endif
    for (int i=0; i<frames; i++) {
        float f1 = dot(bondVector[i], x1[i], 3) * k1;
        float f2 = dot(bondVector[i], x2[i], 3) * k2;
        a1[i]   = f1;
        a2[i]   = f2;
        mid[i]  = (f1+f2)*0.5;
        drop[i] = f2 - f1;
    }
}

void ReadArc::fileHeader(FILE* file, std::string prefix) {
    fprintf(file, " %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s",
    (prefix+"_"+std::to_string(i1)).c_str(),
    (prefix+"avg_"+std::to_string(i1)).c_str(),
    (prefix+"std_"+std::to_string(i1)).c_str(),
    (prefix+"_"+std::to_string(i2)).c_str(),
    (prefix+"avg_"+std::to_string(i2)).c_str(),
    (prefix+"std_"+std::to_string(i2)).c_str(),
    (prefix+"_mid").c_str(),
    (prefix+"avg_mid").c_str(),
    (prefix+"std_mid").c_str(),
    (prefix+"_drop").c_str(),
    (prefix+"avg_drop").c_str(),
    (prefix+"std_drop").c_str()
    );
}

void ReadArc::writeLine(FILE* file, int linen) {
    fprintf(file," %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f",
    a1[linen],
    average(a1,linen),
    std::sqrt(variance(a1,linen)),
    a2[linen],
    average(a2,linen),
    std::sqrt(variance(a2,linen)),
    mid[linen],
    average(mid,linen),
    std::sqrt(variance(mid,linen)),
    drop[linen],
    average(drop,linen),
    std::sqrt(variance(drop,linen))
    );
}

size_t ReadArc::skip2nextAtom(size_t current, std::vector<size_t> &ga1, std::vector<size_t> &ga2, bool smallStep) {
    // Setting a mininum for skipping to reduce cache thrashing
    // but I haven't taken the time to find out how to dynamicall
    // determine the size of the prefetched data
    size_t minSkip = 1024;
    if ((ga1[0] == 0) || (ga2[0] == 0)) return 0;
    size_t d12 = abs(ga1[0]-ga2[0]) - 0;
    if (smallStep and d12 > minSkip) return current+d12;
    if (not smallStep) {
        if (ga1[0] < ga2[0]) { 
            if (ga1[1] == 0) return 0;   
            size_t D12 = ga1[1] - ga2[0];
            if (D12 > minSkip) return current + D12;
            else return 0;
        }
        else {
            if (ga2[1] == 0) return 0;   
            size_t D12 = ga2[1] - ga1[0];
            if (D12 > 1024) return current + D12;
            else return 0;
        }
    }
    return 0;
}
