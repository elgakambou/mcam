#include "BlackScholesModel.hpp"
#include <cmath>

#include <iostream>
using namespace std;

BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect* sigma, PnlVect* divid, PnlVect* spot)
{
    size_ = size;
    r_ = r;
    rho_ = rho;
    sigma_ = sigma;
    divid_ = divid;
    spot_ = spot;

    //Creation de Gamma et de sa factorisation de Cholesky
    L_ =  pnl_mat_create_from_scalar(size_, size_, rho); // A LIBERER 
    pnl_mat_set_diag(L_, 1., 0);
    pnl_mat_chol(L_);
}

void BlackScholesModel::asset(PnlMat *path, double T, int dates, PnlRng *rng)
{
    pnl_mat_set_row(path, spot_, 0);
    PnlVect* norm = pnl_vect_create(path->n);
    PnlVect* L_j = pnl_vect_create(L_->n);
    double exponent_deterministic, exponent_random, sigma_j;
    double time_step = T / (double) dates;
    double sqrt_time_step = sqrt(time_step);
    for (int i=1; i<path->m; i++)
    {
        pnl_vect_rng_normal(norm, norm->size, rng);
        for (int j=0; j<path->n; j++)
        {
            pnl_mat_get_row(L_j, L_, j);
            sigma_j = pnl_vect_get(sigma_, j);
            exponent_deterministic = (r_ - pnl_vect_get(divid_, j) - sigma_j*sigma_j*0.5) * time_step;
            exponent_random =  sqrt_time_step * sigma_j * pnl_vect_scalar_prod(L_j, norm);
            pnl_mat_set(path, i, j, pnl_mat_get(path, i-1, j) * exp(exponent_deterministic + exponent_random));
        }
    }
    pnl_vect_free(&norm);
    pnl_vect_free(&L_j);
}
