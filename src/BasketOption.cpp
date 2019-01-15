#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include <iostream>

#include "BasketOption.hpp"

BasketOption::BasketOption() {
}

BasketOption::BasketOption(double T, double nbTimeSteps, double size, PnlVect *coefs, double strike) {
    coefs_ = pnl_vect_copy(coefs); // taille n
    strike_ = strike; // taille size_
    this->T = T;
    this->nbTimeSteps = nbTimeSteps;
    this->size = size;
    payoffVectMemSpace_ = pnl_vect_create(size);
}

BasketOption::~BasketOption() {
    pnl_vect_free(&coefs_);
    pnl_vect_free(&payoffVectMemSpace_);
}

double BasketOption::payoff(const PnlMat *path) {
    //path de taille [2, size]
    pnl_mat_get_row(payoffVectMemSpace_, path, 1);
    pnl_vect_mult_vect_term(payoffVectMemSpace_, coefs_);

    double d = pnl_vect_sum(payoffVectMemSpace_);
    return (d - strike_ > 0) ? (d - strike_) : 0;
}
