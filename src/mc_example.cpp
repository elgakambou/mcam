#include <iostream>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "PricingResults.hpp"

using namespace std;

int main()
{
    PnlVect *G = pnl_vect_new();
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    int M = 1E5;
    int dim = 2;
    pnl_rng_sseed(rng, time(NULL));

    double acc = 0.;

    for (int i = 0; i < M; i++)
    {
        pnl_vect_rng_normal(G, dim, rng);
        double tmp = pnl_vect_norm_two(G);
        acc += tmp;
    }

    acc /= M;

    cout << PricingResults(acc) << endl;

    pnl_vect_free(&G);
    pnl_rng_free(&rng);
    return 0;
}
