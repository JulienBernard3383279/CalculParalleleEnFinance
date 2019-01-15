#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include <iostream>

#include "PerformanceOption.hpp"

PerformanceOption::PerformanceOption(double T, double nbTimeSteps, double size, PnlVect *coefs, double strike) {
    this->T = T;
    this->nbTimeSteps = nbTimeSteps;
    this->size = size;
    coefs_ = pnl_vect_copy(coefs); // taille size
    payoffVectMemSpaceSums_ = pnl_vect_create(nbTimeSteps+1);
    payoffVectMemSpaceNumerateur_ = pnl_vect_create(nbTimeSteps);
    payoffVectMemSpaceDenominateur_ = pnl_vect_create(nbTimeSteps);
}

double PerformanceOption::payoff(const PnlMat *path) {
    pnl_mat_mult_vect_inplace(payoffVectMemSpaceSums_, path, coefs_);
    pnl_vect_extract_subvect(payoffVectMemSpaceDenominateur_, payoffVectMemSpaceSums_, 0, nbTimeSteps);
    pnl_vect_extract_subvect(payoffVectMemSpaceNumerateur_, payoffVectMemSpaceSums_, 1, nbTimeSteps);
    pnl_vect_div_vect_term(payoffVectMemSpaceNumerateur_,payoffVectMemSpaceDenominateur_);
    pnl_vect_minus_scalar(payoffVectMemSpaceNumerateur_,1.0);

    for (int i=0; i<nbTimeSteps; i++) {
	    if ( pnl_vect_get(payoffVectMemSpaceNumerateur_, i) <= 0.0 ) {
	        pnl_vect_set(payoffVectMemSpaceNumerateur_, i, 0.0);
	    }
    }

    double d = 1 + pnl_vect_sum(payoffVectMemSpaceNumerateur_);
    return d;
}

PerformanceOption::PerformanceOption() {
}

PerformanceOption::~PerformanceOption() {
    pnl_vect_free(&coefs_);
    pnl_vect_free(&payoffVectMemSpaceSums_);
    pnl_vect_free(&payoffVectMemSpaceNumerateur_);
    pnl_vect_free(&payoffVectMemSpaceDenominateur_);
}
