#include "PolynomialRegression.hpp"
#include "iostream"

PolynomialRegression :: PolynomialRegression(int maxDegree,  int sizeOfX) {
    this->id = PNL_BASIS_CANONICAL;
    this->B = pnl_basis_create_from_degree(this->id, maxDegree, sizeOfX); // FAIRE UN DESTRUCTEUR PLUS TARD
    this->maxDegree = maxDegree;
    this->sizeOfX = sizeOfX;  
}


void PolynomialRegression ::  estimate (PnlMat* dataX, PnlVect* dataY, PnlVect* coeff)  {
    pnl_basis_fit_ls (this->B, coeff, dataX, dataY);
}

 double PolynomialRegression :: eval(const PnlVect* coeff, const PnlVect* x) {
            return pnl_basis_eval_vect (this->B, coeff, x) ;
 }