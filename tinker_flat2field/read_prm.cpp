#include "read_prm.hpp"

ReadPrm::ReadPrm(std::string &prmName, std::string &crdName, int &a1, int &a2)
{
    atom = std::vector<int> (2,-1);
    type = std::vector<int> (2,-1);
    alpha = std::vector<float> (2, -100);
    q = std::vector<float> (2, -100);
    atom[0] = a1; atom[1] = a2;
    getTypes(crdName);
    getPrm(prmName);
}

ReadPrm::~ReadPrm(void){}

void ReadPrm::getTypes(std::string &crdName)
{
    int extStart = crdName.rfind(".") + 1;
    int length = crdName.size();
    std::string ext = crdName.substr(extStart,length - extStart);
    if ( strcmp(ext.c_str(),"h5") == 0) {
        const H5std_string FILE_NAME(crdName);
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
        try {
            H5::Group grp = file->openGroup("t_atoms");
            H5::DataSet dataset = grp.openDataSet("type");
            H5::DataSpace dataspace = dataset.getSpace();
            const int one = 1;
            hsize_t dims[1] = { 1 };
            H5::DataSpace mspace2(one, dims);
            hsize_t a1_offset[1] = { (hsize_t)atom[0]-1 }; // offset by -1 b/c .xyz numbering starts
            hsize_t a2_offset[1] = { (hsize_t)atom[1]-1 }; // at 1 while array indexing starts at 0
            hsize_t count[1] = { 1 };
            int t1[1], t2[1];
            dataspace.selectHyperslab(H5S_SELECT_SET, count, a1_offset);
            dataset.read(t1, H5::PredType::NATIVE_INT, mspace2, dataspace);
            dataspace.selectHyperslab(H5S_SELECT_SET, count, a2_offset);
            dataset.read(t2, H5::PredType::NATIVE_INT, mspace2, dataspace);
            type[0] = t1[0];
            type[1] = t2[0];
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
        file->close();
    }
    else {
        int ntype = 0;
        std::string line;
        std::ifstream file(crdName);
        if (file.is_open()) {
            while (file.good() && ntype < 2) {
                std::getline(file, line);
                if (not line.empty()) {
                    std::stringstream linestream(line);
                    std::vector<std::string> lineVector;
                    std::string column;
                    while ( linestream >> column ) {
                        lineVector.push_back(column);
                    }
                    if (lineVector.size() > 5) {
                        if (stod(lineVector[0]) == stoi(lineVector[0])) {
                            if (stoi(lineVector[0]) == atom[0]) {
                                type[0] = stoi(lineVector[5]);
                                ntype++;
                            }
                            if (stoi(lineVector[0]) == atom[1]) {
                                type[1] = stoi(lineVector[5]);
                                ntype++;
                            }
                        }
                    }
                }
            }
        }
        else {
            fprintf(stderr,"\nError: Cannot open %s.\n",crdName.c_str());
            exit(1);
        }
    }
    for (int i=0; i<(int)type.size();i++){
        if (type[i] == -1) {
            fprintf(stderr,"\nError: Cannot find the type for atom %d.\n",atom[i]);
            exit(1);
        }
    }
}

void ReadPrm::getPrm(std::string &prmName)
{
    int nalpha = 0, nq = 0;
    std::string line;
    std::ifstream file(prmName);
    if (file.is_open()) {
        while (file.good() && (nalpha < 2 || nq < 2)) {
            std::getline(file, line);
            if (not line.empty()) {
                std::stringstream linestream(line);
                std::vector<std::string> lineVector;
                std::string column;
                while ( linestream >> column ) {
                    lineVector.push_back(column);
                }
                for (int i=0; i<(int)atom.size(); i++){
                    if (strcmp("multipole", lineVector[0].c_str()) == 0) {
                        if (stoi(lineVector[1]) == type[i]) {
                            q[i] = stod(lineVector[lineVector.size()-1]);
                            nq++;
                        }
                    }
                    if (strcmp("polarize", lineVector[0].c_str()) == 0) {
                        if (stoi(lineVector[1]) == type[i]) {
                            alpha[i] = stod(lineVector[2]);
                            nalpha++;
                        }
                    }
                }
            }
        }
    }
    else {
        fprintf(stderr,"\nError: Cannot open %s.\n",prmName.c_str());
        exit(1);
    }
    for (int i=0; i<(int)type.size();i++){
        if (q[i] == -100) {
            fprintf(stderr,"\nError: Cannot find the charge for atom %d.\n",atom[i]);
            exit(1);
        }
        if (alpha[i] == -100) {
            fprintf(stderr,"\nError: Cannot find the polarizability for atom %d.\n",atom[i]);
            exit(1);
        }

    }
}

int ReadPrm::Atom(int i)
{
    if (i<2) {
        return atom[i];
    }
    else {
        fprintf(stderr,"\nError returning atom index from ReadPrm::Atom.\n");
        exit(1);
    }
}
int ReadPrm::Type(int i)
{
    if (i<2) {
        return type[i];
    }
    else {
        fprintf(stderr,"\nError returning atom type from ReadPrm::Type.\n");
        exit(1);
    }
}
float ReadPrm::Alpha(int i)
{
    if (i<2) {
        return alpha[i];
    }
    else {
        fprintf(stderr,"\nError returning atom polarizability from ReadPrm::Alpha.\n");
        exit(1);
    }
}
float ReadPrm::Q(int i)
{
    if (i<2) {
        return q[i];
    }
    else {
        fprintf(stderr,"\nError returning atom charge from ReadPrm::Q.\n");
        exit(1);
    }
}
