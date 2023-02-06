
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

    //PnlVect* taus_k_1 = pnl_vect_create_from_scalar(nbSamples, T); // les tau_k+1
    PnlVect* taus_k_1 = pnl_vect_create_from_scalar(nbSamples, nbDates); // les tau_k+1

    PnlMat* dataX = pnl_mat_create (nbSamples, mod_->size_); // vecteurs des Stk
    PnlVect* dataY = pnl_vect_create(nbSamples); // vecteurs des payoff
    PnlVect* coeff = pnl_vect_create(estimator->B->nb_func); // vecteur alpha_k 
    PnlVect* S_tk_l = pnl_vect_create(mod_->size_); // les S_tk^l utilisés pour construire dataX
    PnlMat* path_l;
    PnlVect* S_tau_k_l_1 = pnl_vect_create(mod_->size_); 

    double tau_k_l_1;
    double timeStep = T / nbDates;
    double payoff_tk_l;
    for (int k = nbDates - 1; k >=1 ; k--)
    {
        // Construction de dataX, dataY:
        int l = 0; // le nombre de data dans la monnaie
        for (int l = 0; l<nbSamples; l++)
        {
            tau_k_l_1 = pnl_vect_get(taus_k_1, l);
            path_l = (PnlMat*) pnl_array_get(paths, l);
            pnl_mat_get_row(S_tau_k_l_1, path_l, tau_k_l_1); // (int) (tau_k_l_1 / timeStep)
            double payoff = exp(-mod_->r_ * tau_k_l_1*timeStep) * opt_->payoff(S_tau_k_l_1); // tau_k_l_1         
            pnl_mat_get_row(S_tk_l, path_l, k);
            pnl_mat_set_row(dataX, S_tk_l, l);      
            pnl_vect_set(dataY, l, payoff);         
        }
        // //std::cout << compteur << " "  << nbSamples << "\n";
        // PnlMat* dataX_on = pnl_mat_create (compteur, mod_->size_); // vecteurs des Stk
        // PnlVect* dataY_on = pnl_vect_create(compteur); // vecteurs des payoff
        // pnl_mat_extract_subblock (dataX_on, dataX, 0, compteur, 0, mod_->size_);
        // pnl_vect_extract_subvect(dataY_on, dataY, 0, compteur);
       
        
        // Optimisation:
        estimator->estimate(dataX, dataY, coeff);
        //pnl_vect_print(coeff);

        // MAJ des taus_k:
        for (int l=0; l<nbSamples; l++)
        {
            path_l = (PnlMat*) pnl_array_get(paths, l);
            pnl_mat_get_row(S_tk_l, path_l, k);
            payoff_tk_l = exp(mod_->r_* timeStep) *  exp(-mod_->r_ * k * timeStep) * opt_->payoff(S_tk_l); // actualisé
            if (payoff_tk_l >= std::max(estimator->eval(coeff, S_tk_l), 0.0))
                pnl_vect_set(taus_k_1, l, k ); // (double) k * timeStep
        }
    }
    
    double tau1 = 0.;
    double payoff = 0.0;
    for (int l=0; l<nbSamples; l++)
    {
        path_l = (PnlMat*) pnl_array_get(paths, l);
        tau1 = pnl_vect_get(taus_k_1, l); // indice
        pnl_mat_get_row(S_tau_k_l_1, path_l, tau1 ); // (int) (tau1 / timeStep)
        payoff += exp(-mod_->r_ * tau1 * timeStep) * opt_->payoff(S_tau_k_l_1);
        std :: cout << "payoff : " << payoff << "\n";
    }
    pnl_vect_free(&taus_k_1);
    pnl_vect_free(&dataY);
    pnl_vect_free(&coeff);
    pnl_vect_free(&S_tk_l);
    pnl_vect_free(&S_tau_k_l_1);
    pnl_mat_free(&dataX);
    pnl_array_free(&paths);

    double Z0 = opt_->payoff(mod_->spot_);
    double esp  = payoff / nbSamples;
    std :: cout << "Z0 : " << Z0 << " esp : " << esp << "\n";
    return std::max(opt_->payoff(mod_->spot_), payoff / nbSamples);
}



double MonteCarlo :: price2(double T, int nbDates, PnlRng* rng) {

    // les M trajectoires
    PnlArray* paths = pnl_array_create(nbSamples);
    for (int l=0; l<nbSamples; l++)
    {
        PnlMat* path = pnl_mat_create (nbDates + 1, mod_->size_);
        mod_->asset(path, T, nbDates, rng);
        pnl_array_set(paths, l, (PnlObject*) path);
    }

    double step_date = T/ (double) nbDates;
    double Z0 = opt_->payoff(mod_->spot_);
    double tau_k_1 = 0.0;
    PnlVect* tau_tmp_i = pnl_vect_create_from_scalar(nbSamples, T); 
    PnlMat* dataX = pnl_mat_create (nbSamples, mod_->size_); // vecteurs des Stk
    PnlVect* dataY = pnl_vect_create(nbSamples); // vecteurs des payoff
    PnlVect* coeff = pnl_vect_create(estimator->B->nb_func); // vecteur alpha_k 
    PnlVect* S_tk_l = pnl_vect_create(mod_->size_); // les S_tk^l utilisés pour construire dataX
    PnlMat* path_tmp;
    PnlVect* spots_tmp = pnl_vect_create(mod_->size_); 

    for (int k = nbDates - 1; k >=1 ; k--) {

        // Construction de dataX, dataY:
        //pnl_vect_print(tau_tmp_i);
        for (int l = 0; l<nbSamples; l++)
        {
            tau_k_1 = pnl_vect_get(tau_tmp_i, l); // tau_k+1
            path_tmp = (PnlMat*) pnl_array_get(paths, l);
            pnl_mat_get_row(S_tk_l, path_tmp, k);
            pnl_mat_set_row(dataX, S_tk_l, l);
            pnl_mat_get_row(spots_tmp, path_tmp, (int) tau_k_1 / step_date);
            pnl_vect_set(dataY, l, exp(-mod_->r_ * tau_k_1) * opt_->payoff(spots_tmp));
        }

        // evaluation:
        estimator->estimate(dataX, dataY, coeff);

        // MAJ des taus_k:
        for (int l=0; l<nbSamples; l++)
        {
            // evaluation 
            path_tmp = (PnlMat*) pnl_array_get(paths, l);
            pnl_mat_get_row(S_tk_l, path_tmp, k);
            double esp_cond = estimator->eval(coeff, S_tk_l);
            double payoff_capitalized = exp(-mod_->r_ * (double) k * step_date) * opt_->payoff(S_tk_l);
            if ( payoff_capitalized >= std::max(esp_cond, 0.0) ) {
                pnl_vect_set(tau_tmp_i, l, (double) (k*step_date) );
            }
        }
    }
    // calcul du prix 
    double prix = 0.0;
    double tau1 = 0.0;
    pnl_vect_print(tau_tmp_i);
    for (int l=0; l<nbSamples; l++) 
    {
        path_tmp = (PnlMat*) pnl_array_get(paths, l);
        tau1 = pnl_vect_get(tau_tmp_i, l);
        pnl_mat_get_row(spots_tmp, path_tmp, (int) (tau1 / step_date) );
        prix += exp(-mod_->r_ * tau1) * opt_->payoff(spots_tmp);
        //std::cout << prix << "\n";
    }
    // prix = std::max (Z0, prix/ (double) nbSamples);
    // std::cout << Z0 << "\n";
    return std::max (Z0, prix/ (double) nbSamples);;
}
