
#include "MonteCarlo.hpp"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_array.h"


MonteCarlo::MonteCarlo(BlackScholesModel * mod, Option* opt, Quantifier* estimator,  int nbSamples) {
    this->mod_ = mod;
    this->opt_ = opt;
    this->estimator = estimator;
    this->nbSamples = nbSamples;
}

double MonteCarlo :: price(double T, int nbDates, PnlRng* rng)
{
    // Génération des M trajectoires
    PnlArray* paths = pnl_array_create(nbSamples);
    for (int l=0; l<nbSamples; l++)
    {
        PnlMat* path = pnl_mat_create (nbDates + 1, mod_->size_);
        mod_->asset(path, T, nbDates, rng);
        pnl_array_set(paths, l, (PnlObject*) path);
    }

    PnlVect* taus_k_1 = pnl_vect_create_from_scalar(nbSamples, T); // les tau_k+1
    PnlMat* dataX = pnl_mat_create (nbSamples, mod_->size_); // vecteurs des Stk
    PnlVect* dataY = pnl_vect_create(nbSamples); // vecteurs des payoff
    PnlVect* coeff = pnl_vect_create(estimator->B->nb_func); // vecteur alpha_k 
    PnlVect* S_tk_l = pnl_vect_create(mod_->size_); // les S_tk^l utilisés pour construire dataX
    PnlMat* path_l;
    PnlVect* S_tau_k_l_1 = pnl_vect_create(mod_->size_); 

    double tau_k_l_1, timeStep = T / (double) nbDates;
    double payoff_tk_l;
    for (int k = nbDates - 1; k >=1 ; k--)
    {
        // Construction de dataX, dataY:
        for (int l = 0; l<nbSamples; l++)
        {
            tau_k_l_1 = pnl_vect_get(taus_k_1, l);
            path_l = (PnlMat*) pnl_array_get(paths, l);
            pnl_mat_get_row(S_tk_l, path_l, k);
            pnl_mat_set_row(dataX, S_tk_l, l);
            pnl_mat_get_row(S_tau_k_l_1, path_l, (int) (tau_k_l_1 / timeStep));
            pnl_vect_set(dataY, l, exp(-mod_->r_ * tau_k_l_1) * opt_->payoff(S_tau_k_l_1));
        }
        
        // Optimisation:
        estimator->estimate(dataX, dataY, coeff);

        // MAJ des taus_k:
        for (int l=0; l<nbSamples; l++)
        {
            path_l = (PnlMat*) pnl_array_get(paths, l);
            pnl_mat_get_row(S_tk_l, path_l, k);
            payoff_tk_l = exp(-mod_->r_ * (double) k * timeStep) * opt_->payoff(S_tk_l);
            if (payoff_tk_l >= std::max(estimator->eval(coeff, S_tk_l), 0.0))
                pnl_vect_set(taus_k_1, l, (double) k * timeStep);
        }
    }
    
    double tau1, payoff=0.0;
    for (int l=0; l<nbSamples; l++)
    {
        path_l = (PnlMat*) pnl_array_get(paths, l);
        tau1 = pnl_vect_get(taus_k_1, l);
        pnl_mat_get_row(S_tau_k_l_1, path_l, (int) (tau1 / timeStep));
        payoff += exp(-mod_->r_ * tau1) * opt_->payoff(S_tau_k_l_1);
    }

    pnl_vect_free(&taus_k_1);
    pnl_vect_free(&dataY);
    pnl_vect_free(&coeff);
    pnl_vect_free(&S_tk_l);
    pnl_vect_free(&S_tau_k_l_1);
    pnl_mat_free(&dataX);
    pnl_array_free(&paths);

    return std::max(opt_->payoff(mod_->spot_), payoff / (double) nbSamples);
}
