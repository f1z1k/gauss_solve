#include <iostream>
#include <cmath>
#include <cstring>
#include <QDebug>

#include "initu.h"

using namespace std;

// public
TInitU::TInitU(int N, double *x, double *aU, int maxIter) :N(N), maxIter(maxIter), M(maxIter) {
    this->x = new double [N];
    this->aU = new double [N];
    this->prevW = new double [maxIter];
    this->prevA = new double [maxIter];
    this->prevR = new double [maxIter];
    this->curW = new double [maxIter];
    this->curA = new double [maxIter];
    this->curR = new double [maxIter];
    this->gradW = new double [maxIter];
    this->gradA = new double [maxIter];
    this->gradR = new double [maxIter];
    memcpy(this->x, x, sizeof(double) * N);
    memcpy(this->aU, aU, sizeof(double) * N);
    srand(time(NULL));
    findMinE();
}

TInitU::~TInitU() {
    delete[] this->x;
    delete[] this->aU;
    delete[] this->prevW;
    delete[] this->prevA;
    delete[] this->prevR;
    delete[] this->curW;
    delete[] this->curA;
    delete[] this->curR;

    delete[] this->gradW;
    delete[] this->gradA;
    delete[] this->gradR;
}

double TInitU::get(double x) {
    return U(curW, curA, curR, x);
}

double TInitU::findMinE() {
    double minE = 100;
    for (M = 1; M <= maxIter; M++) {
        addWAR();
        double prevMinE = minE;
        minE = findMinEm();
        if (minE >= prevMinE) {
            break;
        }
    }
    return minE;
}

double normalRand() {
    double sum = 0;
    for (int i = 0; i < 12; i++) {
        sum += double(rand()) / RAND_MAX;
    }
    return sum / 12;

    //return rand() / RAND_MAX;
}

double TInitU::T(double c, double prevT) {
    return c * prevT;
}

double TInitU::findMinEm() {
    double curE = E(), prevE;
    deltaE = curE;
    swap(prevW, curW);
    swap(prevA, curA);
    swap(prevR, curR);
    int i = 0;
    double t = T(1, 1);
    while (t > 1e-10) {
        newWAR();
        prevE = curE;
        curE = E();
        deltaE = curE - prevE;
        if (deltaE <= 0 || normalRand() < exp(-deltaE / t)) {
            swap(prevW, curW);
            swap(prevA, curA);
            swap(prevR, curR);
            t = T(0.99, t);
        }
        qDebug() << curE;
    }
    i++;
    t = T(0.99, t);
    return curE;
}

// private

double TInitU::U(double *w, double *a, double *r, double x) {
    double res = 0;
    for (int j = 0; j < M; j++) {
        res += gaussian(w[j], a[j], r[j], x);
    }
    return res;
}

double TInitU::E() {
    double res = 0;
    for (int i = 0; i < N; i++) {
        double d = U(curW, curA, curR, x[i]) - aU[i];
        res += d * d;
    }
    return res;
}

double TInitU::gaussian(double w, double a, double r, double x) {
	double d = x - a;
	return w * exp(-(d * d) / (r * r));
}


void TInitU::addWAR() {
    double a = 0, b = 1;
    curW[M - 1] = double (rand()) / RAND_MAX * (b - a) + a;
    curA[M - 1] = double (rand()) / RAND_MAX * (b - a) + a;
    curR[M - 1] = double (rand()) / RAND_MAX * (b - a) + a;
}

void TInitU::newWAR() {
   for (int l = 0; l < M; l++) {
       gradW[l] = 0;
       gradA[l] = 0;
       gradR[l] = 0;
   }
   double T = 0.00001;
   for (int i = 0; i < N; i++) {
        double delta = U(prevW, prevA, prevR, x[i]) - aU[i];
        for (int l = 0; l < M; l++) {
            double w = prevW[l], a = prevA[l], r = prevR[l];
            double g = gaussian(w, a, r, x[i]);
            double d = x[i] - a;
            gradW[l] += 2 * delta * gaussian(1, a, r, x[i]);
            gradA[l] += 4 * delta * g * d * a / (r * r);
            gradR[l] += 4 * delta * g * d * d / (r * r * r);
        }
    }

    for (int l = 0; l < M; l++) {
        curW[l] = prevW[l] - T * gradW[l];
        curA[l] = prevA[l] - T * gradA[l];
        curR[l] = prevR[l] - T * gradR[l];
    }
}
