/*-------------------------------------------------------------
* This is an abstract Pseudo class
*------------------------------------------------------------*/

#ifndef PSEUDO_H_
#define PSEUDO_H_

#include <iostream>
#include <cassert>
#include <cstdio>
#include "SimulationBox.h"
#include "PolymerChain.h"

class Pseudo
{
protected:
    SimulationBox *sb;
    PolymerChain *pc;

    int n_complex_grid;
    double *expf, *expf_half;

    void set_exp_factor(
        std::array<int,3> nx, std::array<double,3> dx, double ds);
public:
    Pseudo(SimulationBox *sb, PolymerChain *pc);
    virtual ~Pseudo();

    virtual void find_phi(
        double *phi_a,  double *phi_b,
        double *q1_init, double *q2_init,
        double *w_a, double *w_b, double &single_partition) = 0;

    virtual void get_partition(
        double *q1, int n1,
        double *q2, int n2) = 0;

    virtual void update(){
        set_exp_factor(sb->get_nx(), sb->get_dx(), pc->get_ds());
    }

    // Methods for SWIG
    void find_phi(
        double **phi_a, int *len_p_a,
        double **phi_b, int *len_p_b,
        double *q1_init, int len_q1,
        double *q2_init, int len_q2,
        double *w_a, int len_w_a,
        double *w_b, int len_w_b,
        double &single_partition)
    {

        assert(len_q1  == sb->get_n_grid());
        assert(len_q2  == sb->get_n_grid());
        assert(len_w_a == sb->get_n_grid());
        assert(len_w_b == sb->get_n_grid());
        
        double *phi_a_dynamic = (double *) malloc(sb->get_n_grid()*sizeof(double));
        double *phi_b_dynamic = (double *) malloc(sb->get_n_grid()*sizeof(double));
        
        assert(phi_a_dynamic != NULL);
        assert(phi_b_dynamic != NULL);
        
        if (phi_a_dynamic == NULL)
            std::cout << "Failed malloc() for phi_a" << std::endl;
        if (phi_b_dynamic == NULL)
            std::cout << "Failed malloc() for phi_b" << std::endl;
            
        *phi_a = phi_a_dynamic;
        *phi_b = phi_b_dynamic;
        
        *len_p_a = sb->get_n_grid();
        *len_p_b = sb->get_n_grid();
        
        find_phi(*phi_a, *phi_b, q1_init, q2_init, w_a, w_b, single_partition);
    }
    void get_partition(
        double **q1_out, int *len_q1,  
        double **q2_out, int *len_q2, 
        int n1, int n2)
    {
        
        double *q1_out_dynamic = (double *) malloc(sb->get_n_grid()*sizeof(double));
        double *q2_out_dynamic = (double *) malloc(sb->get_n_grid()*sizeof(double));
        
        assert(q1_out_dynamic != NULL);
        assert(q2_out_dynamic != NULL);
        
        if (q1_out_dynamic == NULL)
            std::cout << "Failed malloc() for q1_out" << std::endl;
        if (q2_out_dynamic == NULL)
            std::cout << "Failed malloc() for q2_out" << std::endl;
            
        *q1_out = q1_out_dynamic;
        *q2_out = q2_out_dynamic;
        
        *len_q1 = sb->get_n_grid();
        *len_q2 = sb->get_n_grid();
        
        get_partition(*q1_out, n1, *q2_out, n2);
    }
};
#endif
