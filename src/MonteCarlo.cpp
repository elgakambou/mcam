
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


double MonteCarlo :: price(double T, int nbDates, PnlRng* rng) {

    // data
    int nbFunct = estimator->B->nb_func;
    double stepDate = T/nbDates;
    PnlMat* pathMC = pnl_mat_create (mod_->size_ ,nbDates + 1); // pour la boucle MonteCarlo
    PnlMat* path = pnl_mat_create (mod_->size_ ,nbDates + 1); // pour contenir les trajectoires de S

    //
    double tau1 = 0.0;
    double Un = 0.0;  //
    double U = 0;
    double sommeU = 0;
    double sommek = 0.0;
    double somme = 0.0;


    // boucle MONTE CARLO POUR LE PRIX 
    for (int l1 = 0 ; l1 < nbSamples; l1 ++ ) {
        mod_->asset(pathMC, T, nbDates, rng);
        double Z0 = opt_->payoff(pathMC, 0.0); // a calculer une fois pour optimiser
        double payoff1 = 0.;
 /***************************************************************************************/       
        // MONTE CARLO  sur Tau1 : on calcule plusieurs fois Tau1
        //for (int n = 0; n < nbSamples; n++)  { // nbSample3

            PnlVect* Stk = pnl_vect_create(mod_->size_);
            double tau = T;
            for (int k = nbDates - 1; k >=1 ; k--) {  // a la fin de cette boucle: on a une valeur pour tau1
                PnlVect* dataY = pnl_vect_create(nbSamples); // vecteurs des payoff // nbSample 2
                PnlMat* dataX = pnl_mat_create (nbSamples, mod_->size_); // vecteurs des Stk  // nbSample 2
                PnlVect* coeff = pnl_vect_create(nbFunct);
                // on simule des trajectoires pour l'approxiamtion
                // calcul des coefficient et des fonctions d'estimation pour en deduire TAU1
                for (int l2 = 0; l2 < nbSamples; l2++) {  // nbSample 2
                        mod_->asset(path, T, nbDates, rng);
                        double payoff2 = opt_->payoff(path, tau); // k*stepDate ou tau ?
                        pnl_vect_set(dataY, l2, payoff2);
                        pnl_mat_get_col(Stk, path,  k); // on extrait Stk
                        pnl_mat_set_row(dataX, Stk, l2);
                }
                // on fait notre regression pour estimer l'esperance conditionnelle
                estimator->estimate (dataX, dataY, coeff); 
                //pnl_vect_print(coeff);
                std::cout << "------" << "\n";
                // mise a jour de Tau
                payoff1 = opt_->payoff(pathMC, k*stepDate);
                pnl_mat_get_col(Stk, pathMC, k); // on extrait Stk dans pathMC; k ou  int tk = (int)(t * dates_ / T_);
                double estimation = estimator->eval(coeff, Stk);
                miseAJourTau(tau, k*stepDate, payoff1, estimation);
                //tau = ( payoff1 >= estimation) ? k*stepDate : tau;
                std::cout <<"tau :" << tau <<" payoff : " << payoff1 << "  esp_Cond : " << estimation << " Z0: " << Z0<< "\n";
                // LIBERER LA MEMOIRE 
            }
    /********************************************************************************/
            // Somme des payoffs au temps d'arret  TAU1
            tau1 = tau;
            //std::cout <<"tau :" << tau1 <<" payoff : " << opt_->payoff(pathMC, tau1) << "\n";
            //std::cout << "tau1 : " << tau1 << "\n";
            U =  exp(-mod_->r_ * tau1) * opt_->payoff(pathMC, tau1);
            sommeU += U;
            //sommek += exp(-mod_->r_ * tau1) * opt_->payoff(pathMC, tau1); // ?actalise
        //}
        // Un  = une premiere approxiamtion de UO
        //Un = std::max (Z0, sommek / nbSamples);  // nbSample3
        // std::cout <<"Un :"<< sommek << "\n";
        // return Un;
        //somme += Un; // moyenne des valeurs de U0
        // MEMOIRE 
    }
    // output 
    //std::cout <<"Un :"<< sommek << "\n";
    return sommeU / nbSamples; // nbSample 1
    //return Un;

    //  MEMOIRE A LIBERER
}



double MonteCarlo :: price2(double T, int nbDates, PnlRng* rng)
{
    // Génération des M trajectoires
    PnlArray* paths = pnl_array_create(nbSamples);
    for (int l=0; l<nbSamples; l++)
    {
        PnlMat* path = pnl_mat_create (nbDates + 1, mod_->size_);
        mod_->asset2(path, T, nbDates, rng);
        pnl_array_set(paths, l, (PnlObject*) path);
    }

    PnlVect* taus_k_1 = pnl_vect_create_from_scalar(nbSamples, T); // les tau_k+1
    PnlMat* dataX = pnl_mat_create (nbSamples, mod_->size_); // vecteurs des Stk
    PnlVect* dataY = pnl_vect_create(nbSamples); // vecteurs des payoff
    PnlVect* coeff = pnl_vect_create(estimator->B->nb_func); // vecteur alpha_k 
    PnlMat* taus = pnl_mat_create_from_scalar(nbDates+1, nbSamples, T); // matrice des tau
    PnlVect* S_tk_l = pnl_vect_create(mod_->size_); // les S_tk^l utilisés pour construire dataX
    PnlVect* payoffs_tk = pnl_vect_create(nbSamples); // les payoffs en tk
    PnlVect* estimations_tk = pnl_vect_create(nbSamples); // les alpha_k * H_k(Stk)

    PnlMat* path_l;
    for (int k = nbDates - 1; k >=1 ; k--)
    {
        // Construction de dataX, dataY:
        for (int l = 0; l<nbSamples; l++)
        {
            path_l = (PnlMat*) pnl_array_get(paths, l);
            pnl_mat_get_row(S_tk_l, path_l, k);
            pnl_mat_set_row(dataX, S_tk_l, l);
            pnl_vect_set(dataY, l, opt_->payoff(path_l, pnl_vect_get(taus_k_1, l)));
            pnl_vect_set(payoffs_tk, l, opt_->payoff(path_l, (double) k * T / (double) nbDates));
        }
        
        // Optimisation:
        estimator->estimate(dataX, dataY, coeff);

        // MAJ des taus_k:
        for (int l=0; l<nbSamples; l++)
        {
            path_l = (PnlMat*) pnl_array_get(paths, l);
            pnl_mat_get_row(S_tk_l, path_l, k); // on récup les Stk^l
            pnl_vect_set(estimations_tk, l, estimator->eval(coeff, S_tk_l)); // on calcule alpha_k H_k(Stk^l) 
            if (pnl_vect_get(payoffs_tk, l) > std::max(pnl_vect_get(estimations_tk, l), 0.0))
                pnl_vect_set(taus_k_1, l, (double)k * T / (double) nbDates);
        }
    }
    double tau1, Un_1, prix=0.0;
    for (int l=0; l<nbSamples; l++)
    {
        path_l = (PnlMat*) pnl_array_get(paths, l);
        tau1 = pnl_vect_get(taus_k_1, l);
        Un_1 = exp(-mod_->r_ * tau1) * opt_->payoff(path_l, tau1);
        prix += std::max (opt_->payoff(path_l, 0.0), Un_1);
    }
    return prix/nbSamples;
}



    void MonteCarlo :: miseAJourTau(double& tau, double tk, double payoff, double Ukplus1) {
        if (payoff >= std :: max( Ukplus1, 0.0) ) {
            tau = tk;
        }
    }