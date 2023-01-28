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
    PnlVect* simul_vect = pnl_vect_create(path->n);
    for (int j=0; j<path->m; j++)
    {
        assetVect(simul_vect, j, T, dates, rng);
        pnl_mat_set_row(path, simul_vect, j);
    }
    pnl_vect_free(&simul_vect);
}


// on place les trajectoires par colonne
// path : (nbTimeSteps+1) x d
void BlackScholesModel::asset2(PnlMat *path, double T, int dates, PnlRng *rng)
{
    PnlVect* simul_vect = pnl_vect_create(path->m);
    for (int j=0; j<path->n; j++)
    {
        assetVect(simul_vect, j, T, dates, rng);
        pnl_mat_set_col(path, simul_vect, j);
    }
    pnl_vect_free(&simul_vect);
}




void BlackScholesModel::assetVect(PnlVect *vect, int j, double T, int dates, PnlRng *rng)
{
    double time_step = T / dates;
    double sqrt_time_step = sqrt(time_step);
    double exponent_deterministic, exponent_random;
    PnlVect* norm = pnl_vect_create(L_->n);

    double sigma_j = pnl_vect_get(sigma_, j);
    PnlVect* L_j = pnl_vect_create(L_->n);
    pnl_mat_get_row(L_j, L_, j);
    
    pnl_vect_set(vect, 0, pnl_vect_get(spot_, j)); // ce serait pas j ??
    for (int i=1; i<vect->size; i++)
    {
        pnl_vect_rng_normal(norm, norm->size, rng);
        exponent_deterministic = (r_ - pnl_vect_get(divid_, j) - sigma_j*sigma_j/2) * time_step; // A OPTIMISER 
        exponent_random =  sqrt_time_step * sigma_j * pnl_vect_scalar_prod(L_j, norm);
        pnl_vect_set(vect, i, pnl_vect_get(vect, i-1) * exp(exponent_deterministic + exponent_random)); // LA AUSSI 
    }
    pnl_vect_free(&L_j);
    pnl_vect_free(&norm);
}
