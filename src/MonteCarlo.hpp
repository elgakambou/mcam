#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "Quantifier.hpp"
#include "pnl/pnl_random.h"

class MonteCarlo
{
public:
    BlackScholesModel *mod_; /*! pointeur vers le modèle */
    Option *opt_; /*! pointeur sur l'option */
    Quantifier * estimator;
    int nbSamples;
    MonteCarlo(BlackScholesModel * mod, Option* opt, Quantifier* estimator, int nbSamples);
    /**
     * Calcule le prix de l'option à la date 0
     *
     * @return valeur de l'estimateur Monte Carlo
     */
    double price(double T, int nbDates, PnlRng* rng);

    void miseAJourTau(double& tau, double tk, double payoff, double Ukplus1);

};


