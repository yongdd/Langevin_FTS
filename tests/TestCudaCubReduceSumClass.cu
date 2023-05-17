#include <thrust/reduce.h>
#include "CudaCommon.h"

class TestCudaCub
{

public:

    const int M{32*32*32};

    // variables for cub reduction sum
    size_t temp_storage_bytes;
    double *d_temp_storage;
    double *d_array;
    double *d_array_sum;

    TestCudaCub()
    {
        gpu_error_check(cudaMalloc((void**)&d_array, sizeof(double)*M));
        gpu_error_check(cudaMalloc((void**)&d_array_sum, sizeof(double)*1));

        temp_storage_bytes = 0;
        d_temp_storage = nullptr;

        std::cout << "temp_storage_bytes " << temp_storage_bytes << std::endl; 
        cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_array, d_array_sum, M);
        std::cout << "temp_storage_bytes " << temp_storage_bytes << std::endl;
    };

    ~TestCudaCub()
    {
        cudaFree(d_array);
        cudaFree(d_array_sum);
        cudaFree(d_temp_storage);
    };

};
int main()
{
    try{
        TestCudaCub testcub = TestCudaCub();
        if (testcub.temp_storage_bytes == 0)
            return -1;
        else
            return 0;
    }
    catch(std::exception& exc)
    {
        std::cout << exc.what() << std::endl;
        return -1;
    }
}