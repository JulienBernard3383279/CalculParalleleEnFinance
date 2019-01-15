#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include <iostream>

#include "AsianOption.hpp"

AsianOption::AsianOption() {
}

AsianOption::AsianOption(double T, double nbTimeSteps, double size, PnlVect *coefs, double strike) {
    coefs_ = pnl_vect_copy(coefs); // taille size
    strike_ = strike;
    this->T = T;
    this->nbTimeSteps = nbTimeSteps;
    this->size = size;
    payoffVectMemSpace_ = pnl_vect_create(size);
}

AsianOption::~AsianOption() {
    pnl_vect_free(&coefs_);
    pnl_vect_free(&payoffVectMemSpace_);
}

double AsianOption::payoff(const PnlMat *path) {
    pnl_mat_sum_vect(payoffVectMemSpace_, path, 'r');
    pnl_vect_mult_vect_term(payoffVectMemSpace_, coefs_);
    pnl_vect_mult_scalar(payoffVectMemSpace_, 1.0/((double)nbTimeSteps+1.0) );

    double d = pnl_vect_sum(payoffVectMemSpace_);
    return (d - strike_ > 0) ? (d - strike_) : 0;
}
