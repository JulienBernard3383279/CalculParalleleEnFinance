 /**
 * \file    MonteCarlo.cpp
 * \author  3A IF - Equipe 1
 * \version 1.0
 * \date    13/09/2017
 * \brief   Implémente la structure de la classe MonteCarlo
 *
 * \details
 */

#include <iostream>
#include <cmath>
#include "mpi.h"

#include "MonteCarlo.hpp"

using namespace std;

MonteCarlo::MonteCarlo() {}

MonteCarlo::MonteCarlo(BlackScholesModel *mod, Option *opt, PnlRng *rng, double fdStep, int nbSamples) {
    mod_ = mod;
    opt_ = opt;
    rng_ = rng;
    fdStep_ = fdStep;
    nbSamples_ = nbSamples;
}

MonteCarlo::~MonteCarlo() {}

void MonteCarlo::price_master(double &prix, double &ic) {

    double mySum = 0;
    double mySquaredSum = 0;
    double var = 0;

    double* theirSums;
    double* theirSquaredSums;

    MPI_Status status;
    int mpiSize;
    MPI_Comm_size( MPI_COMM_WORLD, &mpiSize);

    theirSums = new double[mpiSize];
    theirSquaredSums = new double[mpiSize];

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps+1, mod_->size_);

    mod_->initAsset(opt_->nbTimeSteps);
    for (int i = 0; i < nbSamples_-(mpiSize-1)*(nbSamples_/mpiSize); ++i) {
        mod_->postInitAsset(path, opt_->T, opt_->nbTimeSteps, rng_);

        mySum += opt_->payoff(path);
        mySquaredSum += pow(opt_->payoff(path), 2);
    }

    MPI_Gather(&mySum, 1, MPI_DOUBLE, theirSums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&mySquaredSum, 1, MPI_DOUBLE, theirSquaredSums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 1; i < mpiSize ; i++) {

        mySum+=theirSums[i];
        mySquaredSum+=theirSquaredSums[i];
    }

    prix = mySum/nbSamples_ * exp(-mod_->r_ * opt_->T);

    var = exp(-mod_->r_ * opt_->T * 2)
            * (mySquaredSum/nbSamples_ - pow(mySum/nbSamples_, 2));

    ic = 2 * 1.96 * sqrt(var) / sqrt(nbSamples_);

    // Free memory
    pnl_mat_free(&path);
}

void MonteCarlo::price_master_given_precision(double &prix, double &ic, double precision) {
    double mySum = 0;
    double mySquaredSum = 0;
    double var = 0;

    double* theirSums;
    double* theirSquaredSums;

    MPI_Status status;
    int mpiSize;
    MPI_Comm_size( MPI_COMM_WORLD, &mpiSize);

    theirSums = new double[mpiSize];
    theirSquaredSums = new double[mpiSize];

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps+1, mod_->size_);

    mod_->initAsset(opt_->nbTimeSteps);

    int currentNbSamples = 1000;
    int totalNbSamples = currentNbSamples;
    int samplesToSend = 0;

    do {
        for (int i = 0; i < currentNbSamples-(mpiSize-1)*(currentNbSamples/mpiSize); ++i) {
            mod_->postInitAsset(path, opt_->T, opt_->nbTimeSteps, rng_);

            mySum += opt_->payoff(path);
            mySquaredSum += pow(opt_->payoff(path), 2);
        }

        MPI_Gather(&mySum, 1, MPI_DOUBLE, theirSums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&mySquaredSum, 1, MPI_DOUBLE, theirSquaredSums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int i = 1; i < mpiSize ; i++) {

            mySum+=theirSums[i];
            mySquaredSum+=theirSquaredSums[i];
        }

        //cout << "totalNbSamples : " << totalNbSamples << endl;
        prix = mySum/totalNbSamples * exp(-mod_->r_ * opt_->T);

        var = exp(-mod_->r_ * opt_->T * 2)
            * (mySquaredSum/totalNbSamples - pow(mySum/totalNbSamples, 2))
            / totalNbSamples;

        ic = 2 * 1.96 * sqrt(var);

        if (ic > precision) {
            currentNbSamples = (int) (totalNbSamples * ( (ic/precision)*(ic/precision) - 1 )*0.9 + 40);
            totalNbSamples += currentNbSamples;
            samplesToSend = currentNbSamples/4;
        }
        else {
            samplesToSend=0;
        }
        //cout << "Sending "<<samplesToSend<<" samples to compute to other threads.\n";
        MPI_Bcast(&samplesToSend, 1, MPI_INT, 0, MPI_COMM_WORLD);

    } while (ic > precision);

    cout << "Précision atteinte. "<<totalNbSamples<<" samples fûrent nécessaires.\n";

    // Free memory
    pnl_mat_free(&path);
}

void MonteCarlo::price_slave_given_precision() {

    double sum = 0;
    double squaredSum = 0;
    double var = 0;

    MPI_Status status;
    int mpiSize;
    MPI_Comm_size( MPI_COMM_WORLD, &mpiSize);

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps+1, mod_->size_);

    int samplesToReceive = 1000/mpiSize;

    mod_->initAsset(opt_->nbTimeSteps);

    do {
        sum=0;
        squaredSum=0;

        for (int i = 0; i < samplesToReceive; ++i) {
            mod_->postInitAsset(path, opt_->T, opt_->nbTimeSteps, rng_);

            sum += opt_->payoff(path);
            squaredSum = squaredSum + pow(opt_->payoff(path), 2);
        }

        MPI_Gather(&sum, 1, MPI_DOUBLE, 0, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&squaredSum, 1, MPI_DOUBLE, 0, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Bcast(&samplesToReceive, 1, MPI_INT, 0, MPI_COMM_WORLD);

    } while (samplesToReceive>0);

    // Free memory
    pnl_mat_free(&path);
}



void MonteCarlo::price_slave() {

    double sum = 0;
    double squaredSum = 0;
    double var = 0;

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps+1, mod_->size_);

    mod_->initAsset(opt_->nbTimeSteps);

    for (int i = 0; i < nbSamples_-3*(nbSamples_/4); ++i) {
        mod_->postInitAsset(path, opt_->T, opt_->nbTimeSteps, rng_);

        sum += opt_->payoff(path);
        squaredSum = squaredSum + pow(opt_->payoff(path), 2);
    }

    MPI_Gather(&sum, 1, MPI_DOUBLE, 0, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&squaredSum, 1, MPI_DOUBLE, 0, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Free memory
    pnl_mat_free(&path);
}
