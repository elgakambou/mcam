#include <iostream>
#include <ctime>
#include <string>
#include "jlparser/parser.hpp"
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "BlackScholesModel.hpp"
#include "PerformanceOption.hpp"
#include "BasketOption.hpp"
#include "GeometricOption.hpp"
#include "MonteCarlo.hpp"
#include "Quantifier.hpp"
#include "PolynomialRegression.hpp"
#include "PricingResults.hpp"

using namespace std;

int main(int argc, char **argv)
{
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    // pnl_rng_sseed(rng, 0.0);
    int size, nbTimeSteps = 10;
    int n_samples, maxDegree;
    double r, rho, T, K, strike;
    PnlVect *sigma, *divid, *spot, *lambda;
    string type ;

    // Param *P = new Parser("/home/hanriotr/Scholar/3A/OA/mcam/dat/geom_d2.txt");
    Param *P = new Parser(argv[1]);
    P->extract("maturity", T);
    P->extract("model size", size);
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);
    P->extract("interest rate", r);
    if (P->extract("dividend rate", divid, size, true) == false)
    {
        divid = pnl_vect_create_from_zero(size);
    }
    P->extract("strike", strike);
    P->extract("option type", type);
    P->extract("dates", nbTimeSteps);
    P->extract("MC iterations", n_samples);
    P->extract("degree for polynomial regression", maxDegree);

    // test print 
    // cout << "option type " << type << endl;
    // cout << "maturity " << T << endl;
    // cout << "strike " << strike << endl
    // cout << "model size " << size << endl;
    // cout << "interest rate " << r << endl;
    // cout << "dividend rate ";.
    // pnl_vect_print_asrow(divid);
    // cout << "spot ";
    // pnl_vect_print_asrow(spot);
    // cout << "volatility "
    // pnl_vect_print_asrow(sigma);
    // cout << "MC iterations " << n_samples << endl;
    //PnlMat *path = pnl_mat_create(size, nbTimeSteps+1);
    BlackScholesModel model(size, r, rho, sigma, divid, spot);
    //model.asset(path, T, nbTimeSteps, rng);
    //pnl_mat_print(path);

    GeometricOption geom(T, nbTimeSteps, size, strike);
    // MC
    PolynomialRegression estimator(maxDegree, size); //size
    // pnl_basis_print(estimator.B);
    // exit(0);
    MonteCarlo mc(&model, &geom, &estimator, n_samples); //n_samples

    // l'option 
    if (type == "exchange") {
        P->extract("payoff coefficients", lambda, size);
        BasketOption basket(T, nbTimeSteps, size, strike, lambda);
        mc.opt_ = &basket;
    }
    else if (type == "geometric_put"){
        mc.opt_ = &geom;
    }
    else if (type == "bestof") {
        P->extract("payoff coefficients", lambda, size);
        PerformanceOption perf(T, nbTimeSteps, size, strike, lambda);
        mc.opt_ = &perf;

    }
    else {
        std::cout << " option inconnue" << "\n";
        exit(0);
    }
    double myprice = mc.price(T, nbTimeSteps, rng);
    std::cout << PricingResults(myprice) << std::endl;


    //pnl_mat_free(&path);
    pnl_rng_free(&rng);
    return 0;
}
