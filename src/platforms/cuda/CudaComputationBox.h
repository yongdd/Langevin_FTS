/*-------------------------------------------------------------
* This is a derived CudaComputationBox class
*------------------------------------------------------------*/

#ifndef CUDA_SIMULATION_BOX_H_
#define CUDA_SIMULATION_BOX_H_

#include <vector>
#include "ComputationBox.h"

class CudaComputationBox : public ComputationBox
{
private:
    double *sum, *d_sum;   // temporal storage for reduction in integral_gpu
    double *d_multiple;    // temporal storage for mutiple_inner_product_gpu
    double *d_dv; // dV for GPU

    // variables for cub reduction sum
    size_t temp_storage_bytes = 0;
    double *d_temp_storage = NULL;
    double *d_sum_out;

    void initialize();
public:
    CudaComputationBox(std::vector<int> nx, std::vector<double> lx);
    ~CudaComputationBox() override;

    double integral_gpu(double *d_g);
    double inner_product_gpu(double *d_g, double *d_h);
    double inner_product_inverse_weight_gpu(double *d_g, double *d_h, double *d_w);
    double mutiple_inner_product_gpu(int n_comp, double *d_g, double *d_h);
    void set_lx(std::vector<double> new_lx) override;
};
#endif
