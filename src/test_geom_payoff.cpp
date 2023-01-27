#include <iostream>
#include <ctime>
#include <string>
#include "jlparser/parser.hpp"
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "BlackScholesModel.hpp"
#include "GeometricOption.hpp"


using namespace std;

int main(int argc, char **argv)
{
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    // pnl_rng_sseed(rng, 0.0);
    int size, nbTimeSteps = 10;
    double r, rho, T, K, strike;
    PnlVect *sigma, *divid, *spot;

    // Param *P = new Parser("../dat/geom_d2.txt");
    Param *P = new Parser("../dat/geom_d10.txt");
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

    PnlMat *path = pnl_mat_create(size, nbTimeSteps+1);
    BlackScholesModel model(size, r, rho, sigma, divid, spot);
    model.asset(path, T, nbTimeSteps, rng);
    pnl_mat_print(path);
    
    GeometricOption geom(T, nbTimeSteps, size, strike);

    double t = 0.5;
    double payoff_ = geom.payoff(path, t);

    cout << "GeometricPayoff(" << t << ") = " << payoff_ << endl;

    pnl_mat_free(&path);
    pnl_rng_free(&rng);
    return 0;
}
