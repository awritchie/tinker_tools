#ifndef tinker2field_hpp
#define tinker2field_hpp
#include "read_arc.hpp"
#include "read_prm.hpp"
#include <numeric>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

template <typename T>
float dot(const std::vector<T> &a, const std::vector<T> b, int size) {
    float result = 0;
    for (int i=0; i<size; i++){
        result += a[i] * b[i];
    }
    return result;
}

template <typename T>
float average(const std::vector<T> &a, int size) {
    float sum = 0;
    for (int i=0; i<size+1; i++){
        sum += a[i];
    }
    return sum / (size+1);
}

template <typename T>
float square_average(const std::vector<T> &a, int size) {
    float sum = 0;
    for (int i=0; i<size+1; i++){
        sum += a[i]*a[i];
    }
    return sum / (size+1);
}

template <typename T>
float variance(const std::vector<T> &a, int size) {
    return std::max(0.,square_average(a,size) - std::pow(average(a,size),2));
}
// var = <x^2> - <x>^2
#endif
