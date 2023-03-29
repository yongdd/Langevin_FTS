#define THRUST_IGNORE_DEPRECATED_CPP_DIALECT
#define CUB_IGNORE_DEPRECATED_CPP_DIALECT

#include <iostream>
#include <cstdlib>
#include <string>

#include "CudaCommon.h"

void throw_on_cuda_error(cudaError_t code, const char *file, int line, const char *func)
{
    if (code != cudaSuccess){
        std::string file_and_line("File: \"" + std::string(file) + "\", line: " + std::to_string(line) + ", function <" + std::string(func) + ">");
        throw thrust::system_error(code, thrust::cuda_category(), file_and_line);
    }
}

CudaCommon::CudaCommon()
{
    try{
        // intialize NUM_BLOCKS and NUM_THREADS
        const char *ENV_N_BLOCKS  = getenv("LFTS_GPU_NUM_BLOCKS");
        const char *ENV_N_THREADS = getenv("LFTS_GPU_NUM_THREADS");

        std::string env_var_n_blocks (ENV_N_BLOCKS  ? ENV_N_BLOCKS  : "");
        std::string env_var_n_threads(ENV_N_THREADS ? ENV_N_THREADS : "");

        if (env_var_n_blocks.empty())
            this->n_blocks = 256;
        else
            this->n_blocks = std::stoi(env_var_n_blocks);

        if (env_var_n_threads.empty())
            this->n_threads = 256;
        else
            this->n_threads = std::stoi(env_var_n_threads);

        // the number of GPUs
        int devices_count;
        gpu_error_check(cudaGetDeviceCount(&devices_count));
        const char *ENV_N_GPUS = getenv("LFTS_NUM_GPUS");
        std::string env_var_n_gpus (ENV_N_GPUS  ? ENV_N_GPUS  : "");

        if (env_var_n_gpus.empty())
            n_gpus = 1;
        else
            n_gpus = std::min(std::min(std::stoi(env_var_n_gpus), devices_count), MAX_GPUS);

        // check if can access peer GPUs
        if (n_gpus > 1)
        {
            int can_access_from_0_to_1;
            int can_access_from_1_to_0;
            gpu_error_check(cudaDeviceCanAccessPeer(&can_access_from_0_to_1, 0, 1));
            gpu_error_check(cudaDeviceCanAccessPeer(&can_access_from_1_to_0, 1, 0));

            if (can_access_from_0_to_1 == 1 && can_access_from_1_to_0 == 1)
            {
                gpu_error_check(cudaSetDevice(0));
                gpu_error_check(cudaDeviceEnablePeerAccess(1, 0));
                gpu_error_check(cudaSetDevice(1));
                gpu_error_check(cudaDeviceEnablePeerAccess(0, 0));
            }
            else
            {
                std::cout << "Could not establish peer access between GPUs." << std::endl;
                std::cout << "Only one GPU will be used." << std::endl;
                n_gpus = 1;
            }
        }
        gpu_error_check(cudaSetDevice(0));
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void CudaCommon::set(int n_blocks, int n_threads, int process_idx)
{
    int devices_count;

    this->set_n_blocks(n_blocks);
    this->set_n_threads(n_threads);

    // change GPU setting
    gpu_error_check(cudaGetDeviceCount(&devices_count));
    gpu_error_check(cudaSetDevice(process_idx%devices_count));
}
int CudaCommon::get_n_blocks()
{
    return n_blocks;
}
int CudaCommon::get_n_threads()
{
    return n_threads;
}
int CudaCommon::get_n_gpus()
{
    return n_gpus;
}
void CudaCommon::set_n_blocks(int n_blocks)
{
    this->n_blocks = n_blocks;
}
void CudaCommon::set_n_threads(int n_threads)
{
    this->n_threads = n_threads;
}
void CudaCommon::set_idx(int process_idx)
{
    int devices_count;

    // change GPU setting
    gpu_error_check(cudaGetDeviceCount(&devices_count));
    gpu_error_check(cudaSetDevice(process_idx%devices_count));
}
__global__ void multi_real(double* dst,
                          double* src1,
                          double* src2,
                          double  a, const int M)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < M)
    {
        dst[i] = a * src1[i] * src2[i];
        i += blockDim.x * gridDim.x;
    }
}

__global__ void mutiple_multi_real(int n_comp,
                          double* dst,
                          double* src1,
                          double* src2,
                          double  a, const int M)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < M)
    {  
        dst[i] = a * src1[i] * src2[i];
        for(int n = 1; n < n_comp; n++)
            dst[i] += a * src1[i+n*M] * src2[i+n*M];
        i += blockDim.x * gridDim.x;
    }
}

__global__ void divide_real(double* dst,
                          double* src1,
                          double* src2,
                          double  a, const int M)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < M)
    {
        dst[i] = a * src1[i]/src2[i];
        i += blockDim.x * gridDim.x;
    }
}
__global__ void add_multi_real(double* dst,
                             double* src1,
                             double* src2,
                             double  a, const int M)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < M)
    {
        dst[i] = dst[i] + a * src1[i] * src2[i];
        i += blockDim.x * gridDim.x;
    }
}

__global__ void lin_comb(double* dst,
                        double a,
                        double* src1,
                        double b,
                        double* src2,
                        const int M)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < M)
    {
        dst[i] = a*src1[i] + b*src2[i];
        i += blockDim.x * gridDim.x;
    }
}

__global__ void add_lin_comb(double* dst,
                           double a,
                           double* src1,
                           double b,
                           double* src2,
                           const int M)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < M)
    {
        dst[i] = dst[i] + a*src1[i] + b*src2[i];
        i += blockDim.x * gridDim.x;
    }
}

__global__ void multi_complex_real(ftsComplex* dst,
                                 double* src, const int M)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < M)
    {
        dst[i].x = dst[i].x * src[i];
        dst[i].y = dst[i].y * src[i];
        i += blockDim.x * gridDim.x;
    }
}

__global__ void multi_complex_real(ftsComplex* dst,
                                 double* src, double a, const int M)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < M)
    {
        dst[i].x = a * dst[i].x * src[i];
        dst[i].y = a * dst[i].y * src[i];
        i += blockDim.x * gridDim.x;
    }
}

__global__ void multi_complex_conjugate(double* dst,
                                 ftsComplex* src1,
                                 ftsComplex* src2, const int M)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < M)
    {
        dst[i] = src1[i].x * src2[i].x + src1[i].y * src2[i].y;
        i += blockDim.x * gridDim.x;
    }
}
