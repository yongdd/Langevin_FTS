#ifndef CPU_ARRAY_H_
#define CPU_ARRAY_H_

#include <vector>

#include "Exception.h"
#include "Array.h"

class CpuArray : public Array
{

private:
    double *data = nullptr;
    unsigned int size = 0;
public:
    CpuArray(unsigned new_size);
    CpuArray(double* new_data, unsigned new_size);
    CpuArray(const Array& array);

    ~CpuArray();

    // overloading for array operation
    void add(const Array& src_1, const Array& src_2) override;
    void subtract(const Array& src_1, const Array& src_2) override;
    void multiply(const Array& src_1, const Array& src_2) override;
    void divide(const Array& src_1, const Array& src_3) override;

    // arithmetic operations with a float number
    void linear_scaling(const Array& src, const double a, const double b) override;

    // access element of array
    void operator=(const Array&) override;
    void set_data(double *, unsigned int) override;

    // return array as vector
    std::vector<double> to_vector() const override;

    // return data pointer
    double* get_ptr() const override;

    // return size of data
    unsigned int get_size() const override;
    
    // access element
    double operator[](unsigned int) const override;
};
#endif