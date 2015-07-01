#include "tinkerxyz.hpp"

TinkerXYZ::TinkerXYZ(const std::string &file) : xyz(_xyz), atom(_atom), atype(_atype), bonded(_bonded), box(_box)
{
    filename = file;
    open_tinker_read();
    read_tinker_structure();
    numframes = 0;
    /* For output purposes only  */
    threshold = std::vector<int> {1, 50, 250, 1000, std::numeric_limits<int>::max()};
    increment = std::vector<int> {1, 10,  50,  100}; 
}

TinkerXYZ::~TinkerXYZ(void)
{
    close_tinker_read();    
}

bool TinkerXYZ::isBox(char *c)
{
    int i=0;
    while (c[i]) {
        if (isalpha(c[i])) {
            return false;
        }
        i++;
    }
    return true;
}

void TinkerXYZ::open_tinker_read()
{
    file = fopen(filename.c_str(), "rb");
    if (!file) {
        fprintf(stderr, "TinkerXYZ::open_tinker_read failed to open file '%s'.\n", filename.c_str());
        exit(1); 
    }

    /* First line is the number of atoms */
    int i = fscanf(file, "%d", &numatoms);
    if (i<1) {
        fprintf(stderr, "\nError: Tinker file '%s' should have the number of atoms in the first line.\n", filename.c_str());
        exit(1);
    }    
    /* skip to the end of the line */
    while(getc(file) != '\n');

    /* Determine if simulation box values are present */
    octahedron = false;
    triclinic = false;
    monoclinic = false;
    orthogonal = false;
    char cl[1024], fbuffer[1024], *k;
    float xbox, ybox, zbox, alpha, beta, gamma;
    k = fgets(fbuffer, 1024, file);
    int j = sscanf(fbuffer, "%d %s", &xbox, cl);
    if (isBox(cl)) {
        sscanf(fbuffer, "%f %f %f %f %f %f", &xbox, &ybox, &zbox, &alpha, &beta, &gamma);
        hasBox = true;
        if (alpha == 90 && beta == 90 && gamma == 90) { 
            orthogonal = true;
        } else if (alpha == 90 && gamma == 90) { 
            monoclinic = true;
        } else { 
            triclinic = true;
        }
    } else {
        hasBox = false;
    }

    fprintf(stderr, "TinkerXYZ::open_tinker_read succeeded opening file '%s'.\n", filename.c_str());
    fprintf(stderr, "Number of atoms: %d\n", numatoms);
}

void TinkerXYZ::read_tinker_structure()
{
    int ndx, type, bond;
    std::string atomname;
    float coord;
    char fbuffer[1024], *k;
    
    _atom.reserve(numatoms);
    _atype.reserve(numatoms);
    _bonded.reserve(numatoms);

    for (int i=0; i<numatoms; i++) { 
        k = fgets(fbuffer, 1024, file);
        if (k == NULL) {
            fprintf(stderr, "tinker structure missing atom(s) in file '%s'.\n", filename.c_str());
            fprintf(stderr, "tinker structure expecting '%d' atoms, found only '%d'.\n", numatoms, i+1);
            exit(1);
        } else {
            std::stringstream linestream(fbuffer);
            linestream >> ndx >> atomname >> coord >> coord >> coord >> type;

            std::vector<int> bonds;
            while ( linestream >> bond ) {
                bonds.push_back(bond);
            }
            _atom.push_back(atomname);
            _atype.push_back(type);
            _bonded.push_back(bonds);
        }
    }

    rewind(file);
}

void TinkerXYZ::close_tinker_read()
{
    fclose(file);
}

bool TinkerXYZ::get_next_frame() 
{
    _xyz = std::vector<std::vector<float> > (numatoms, std::vector<float> (3,0));
    char atom_name[1025], fbuffer[1024], *k;
    int j, ndx;
    float x, y, z;
    /* Remove first line, which has the number of atoms */
    if (NULL == fgets(fbuffer, 1024, file)) return false;
    /* If box dimension present, remove second line */
    if (hasBox) {
        if (NULL == fgets(fbuffer, 1024, file)) return false;
        float l1, l2, l3, a1, a2, a3;
        j = sscanf(fbuffer, "%f %f %f %f %f %f", &l1, &l2, &l3, &a1, &a2, &a3);
        if (j != 6) { 
            fprintf(stderr, "unable to read box dimensions.\n");
            fprintf(stderr, "%s", fbuffer);
            return false;
        } else {
            std::vector<float> b = { l1, l2, l3, a1, a2, a3 };
            _box = b;
        }
    }

    /* Read the coordinates */
    for (int i=0; i<numatoms; i++) { 
        k = fgets(fbuffer, 1024, file);
        j = sscanf(fbuffer, "%d %s %f %f %f", &ndx, atom_name, &x, &y, &z);
        if (k == NULL) {
            return false;
        } else if (j < 5) {
            fprintf(stderr, "tinker timestep missing type or coordinate(s) in file '%s' for atom '%d'.\n", filename.c_str(), i+1);
            fprintf(stderr,"%s", fbuffer);
        } else if ( j >= 5) {
            _xyz[i][0] = x;
            _xyz[i][1] = y;
            _xyz[i][2] = z;
        } else {
            break;
        }
    }
    numframes++;
    for (int i=1; i<(int) threshold.size(); i++){
        if ((numframes <= threshold[i]) 
            && (numframes % increment[i-1] ==0)) {
           std::cout << "\rDone reading frame " << numframes << std::flush; // << std::endl; 
           break;
        }
    }

    return true;
}

int TinkerXYZ::natoms() {
    return numatoms;
}

int TinkerXYZ::nframes() {
    return numframes;
}

void TinkerXYZ::isOctahedron() { 
    if (!orthogonal) {
        fprintf(stderr, "WARNING!! Trying to assign 'octahedron' PBC geometry, but the box geometry (alpha, beta, gamma) was not found to be (90.0, 90.0, 90.0).\n");
    }
    octahedron = true;
    orthogonal = false;
}

std::vector<float> TinkerXYZ::pair_minimum_vector(const int &a, const int &b){
    return pair_minimum_vector(_xyz[a], _xyz[b]);
}

std::vector<float> TinkerXYZ::pair_minimum_vector(const std::vector<float> &a, const std::vector<float> &b){
    /*
         Calculates vector between points, taking periodicity 
         into account
    */
    std::vector<float> res (3,0);
    for (int i=0; i<3; i++) { 
        res[i] = b[i] - a[i];
    }
    if (orthogonal) {
        for (int i=0; i<3; i++) {
            while (std::abs(res[i]) > _box[i]*0.5) {
                res[i] -= _box[i] * ((res[i] < 0) ? -1 : (res[i] > 0));
            }
        }
    } else if (monoclinic) {
        float beta_sin = std::sin(_box[4]/RADIAN);
        float beta_cos = std::cos(_box[4]/RADIAN);
        res[2] = res[2] / beta_sin;
        res[0] = res[0] - beta_cos;
        for (int i=0; i<3; i++) {
            while (std::abs(res[i]) > _box[i]*0.5) {
                res[i] -= _box[i] * ((res[i] < 0) ? -1 : (res[i] > 0));
            }
        }
        res[0] = res[0] + res[2] * beta_cos;
        res[2] = res[2] * beta_sin;
    } else if (triclinic) {
        float alpha_cos = std::cos(_box[3]/RADIAN);
        float beta_sin = std::sin(_box[4]/RADIAN);
        float beta_cos = std::cos(_box[4]/RADIAN);
        float gamma_sin = std::sin(_box[5]/RADIAN);
        float gamma_cos = std::cos(_box[5]/RADIAN);
        float beta_term = (alpha_cos - beta_cos * gamma_cos) / gamma_sin;
        float gamma_term = std::sqrt(beta_sin*beta_sin - beta_term*beta_term);
        res[2] = res[2] / gamma_term;
        res[1] = (res[1] - res[2]*beta_term)/ gamma_sin;
        res[0] = res[0] - res[1] * gamma_cos - res[2] * beta_cos;
        for (int i=0; i<3; i++) {
            while (std::abs(res[i]) > _box[i]*0.5) {
                res[i] -= _box[i] * ((res[i] < 0) ? -1 : (res[i] > 0));
            }
        }
        res[0] = res[0] + res[1] * gamma_cos + res[2] * beta_cos;
        res[1] = res[1] * gamma_sin + res[2] * beta_term;
        res[2] = res[2] * gamma_term;
    } else if (octahedron) {
        for (int i=0; i<3; i++) {
            while (std::abs(res[i]) > _box[i]*0.5) {
                res[i] -= _box[i] * ((res[i] < 0) ? -1 : (res[i] > 0));
            }
        }
        if (std::abs(res[0]) + std::abs(res[1]) + std::abs(res[2]) > 0.75 * _box[0]) {
            for (int i=0; i<3; i++) {
                res[i] -= _box[i] * 0.5 * ((res[i] < 0) ? -1 : (res[i] > 0));
            }
        }
    }    
    return res;
}

float TinkerXYZ::distance(const int &a, const int &b) {
    return distance(_xyz[a], _xyz[b]);
}

float TinkerXYZ::distance(const std::vector<float> &a, const std::vector<float> &b) {
    /*
         Calculates distance between points, taking periodicity 
         into account
    */
    std::vector<float> ri = pair_minimum_vector(a, b);
    return std::sqrt(vdot(ri,ri));
}

float TinkerXYZ::angle(const int &a, const int &b, const int &c) {
    return angle(_xyz[a], _xyz[b], _xyz[c]);
}

float TinkerXYZ::angle(const std::vector<float> &a, const std::vector<float> &b, const std::vector<float> c) {
    /*
         Calculates angle between points, taking periodicity 
         into account
    */
    std::vector<float> ba = pair_minimum_vector(b, a);
    std::vector<float> bc = pair_minimum_vector(b, c);
    float rba = std::sqrt(vdot(ba,ba));
    float rbc = std::sqrt(vdot(bc,bc));
    for (int i=0; i<3; i++) {
        ba[i] /= rba;
        bc[i] /= rbc;
    }
    float cos = vdot(ba,bc);
    return std::acos(cos) * RADIAN;
}

float TinkerXYZ::dihedral(const int &a, const int &b, const int &c, const int &d) {
    return dihedral(_xyz[a], _xyz[b], _xyz[c], _xyz[d]);
}

float TinkerXYZ::dihedral(const std::vector<float> &a, const std::vector<float> &b, const std::vector<float> &c, const std::vector<float> &d) {
    /*
         Calculates dihedral angle between points, taking periodicity 
         into account
    */
    std::vector<float> b1 = pair_minimum_vector(a, b);
    std::vector<float> b2 = pair_minimum_vector(b, c);
    std::vector<float> b3 = pair_minimum_vector(c, d);
    float rb1 = std::sqrt(vdot(b1,b1));
    float rb2 = std::sqrt(vdot(b2,b2));
    float rb3 = std::sqrt(vdot(b3,b3));
    for (int i=0; i<3; i++) {
        b1[i] /= rb1;
        b2[i] /= rb2;
        b3[i] /= rb3;
    }
    std::vector<float> n1 = vcross(b1, b2);
    std::vector<float> n2 = vcross(b2, b3);
    std::vector<float> m1 = vcross(n1, b2);
    float x = vdot(n1,n2);
    float y = vdot(m1,n2);
    return - std::atan2(y,x) * RADIAN;
}

float TinkerXYZ::vdot(const std::vector<float> &a, const std::vector<float> &b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


std::vector<float> TinkerXYZ::vcross(const std::vector<float> &a, const std::vector<float> &b) {
    std::vector<float> res (3,0);
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
    return res;
}

float TinkerXYZ::volume() {
    if (orthogonal) { 
        return _box[0]*_box[1]*_box[2];
    } else if (octahedron) { 
        return 0.5*(_box[0]*_box[1]*_box[2]);
    } else if (monoclinic || triclinic) {
        float alpha_cos = std::cos(_box[3]/RADIAN);
        float beta_sin = std::sin(_box[4]/RADIAN);
        float beta_cos = std::cos(_box[4]/RADIAN);
        float gamma_sin = std::sin(_box[5]/RADIAN); 
        float gamma_cos = std::cos(_box[5]/RADIAN);
        float beta_term = (alpha_cos - beta_cos * gamma_cos) / gamma_sin;
        float gamma_term = std::sqrt(beta_sin*beta_sin - beta_term*beta_term);
        return (gamma_sin*gamma_term)*(_box[0]*_box[1]*_box[2]);
    } else {
        fprintf(stderr, "ERROR: Box dimensions undefined, cannot determine box volume.\n");
        return -1;
    }
}

