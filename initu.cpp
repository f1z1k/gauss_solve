#include <iostream>
#include <cmath>
#include <cstring>
#include <QDebug>

#include "initu.h"

using namespace std;

// public
TInitU::TInitU(int N, double *x, double *aU, int maxIter, double eps)
    :N(N), maxIter(maxIter), M(maxIter), eps(eps)
{
    this->x = new double [N];
    this->aU = new double [N];
    this->newW = new double [maxIter];
    this->newA = new double [maxIter];
    this->newR = new double [maxIter];
    this->curW = new double [maxIter];
    this->curA = new double [maxIter];
    this->curR = new double [maxIter];
    this->minW = new double [maxIter];
    this->minA = new double [maxIter];
    this->minR = new double [maxIter];
    memcpy(this->x, x, sizeof(double) * N);
    memcpy(this->aU, aU, sizeof(double) * N);
    srand(time(NULL));
    
    minX = x[0];
    maxX = x[0];
    minU = aU[0];
    maxU = aU[0];
    for (int i = 1; i < N; i++) {
        if (x[i] < minX) {
            minX = x[i];
        } else if (x[i] > maxX) {
            maxX = x[i];
        }
        if (aU[i] < minU) {
            minU = aU[i];
        } else if (aU[i] > maxU) {
            maxU = aU[i];
        }
    }
    findMinE();
}

TInitU::~TInitU() {
    delete[] this->x;
    delete[] this->aU;
    delete[] this->minW;
    delete[] this->minA;
    delete[] this->minR;
    delete[] this->curW;
    delete[] this->curA;
    delete[] this->curR;

    delete[] this->newW;
    delete[] this->newA;
    delete[] this->newR;
}

double TInitU::get(double x) {
    return U(minW, minA, minR, x);
}

// private

double TInitU::findMinE() {
    M = 1;
    firstGaussian();
    minE = E();
    minE = findMinEm();
    qDebug() << minE;

    for (M = 2; M <= maxIter; M++) {
        addGaussian();
        minE = findMinEm();
        qDebug() << minE;
    }

    for (int l = 0; l < maxIter; l++) {
        qDebug() << minW[l] << " " << minA[l] << " " << minR[l];
    }
    return minE;
}

double TInitU::findMinEm() {
    double curE = E();
    t = T(1, 1000);
    while (t > eps) {
        t = T(0.99, t);
        if (curE < minE) {
            minE = curE;
            memcpy(minW, curW, sizeof(double) * M);
            memcpy(minA, curA, sizeof(double) * M);
            memcpy(minR, curR, sizeof(double) * M);
        }
        while (1) {
            g();
            double newE = E();
            double alfa = double(rand()) / RAND_MAX;
            if (alfa < h(newE - curE, t)) {
                swap(curW, newW);
                swap(curA, newA);
                swap(curR, newR);
                curE = newE;
                break;
            }
        }
    }
    return minE;
}

void TInitU::firstGaussian() {
    double t = double (rand()) / RAND_MAX * (maxU - minU) + minU;
    minW[M - 1] = curW[M - 1] = newW[M - 1] = min(max(t, 0.001), maxU);
    t = double (rand()) / RAND_MAX * (maxX - minX) + minX;
    minA[M - 1] = curA[M - 1] = newA[M - 1] = min(max(t, minX), maxX);
    t = double (rand()) / RAND_MAX * (maxX - minX) + minX;
    minR[M - 1] = curR[M - 1] = newR[M - 1] = min(max(t, 0.001), maxX);
}

void TInitU::addGaussian() {
    minW[M - 1] = curW[M - 1] = newW[M - 1] = minU;
    minA[M - 1] = curA[M - 1] = newA[M - 1] = 0.5 * (maxX + minX);
    minR[M - 1] = curR[M - 1] = newR[M - 1] = maxU;
}

double TInitU::T(double c, double prevT) {
    return c * prevT;
}

void TInitU::g() {
    for (int l = 0; l < M; l++) {
        double sumW = 0, sumA = 0, sumR = 0;
        for (int i = 0; i < 12; i++) {
            sumW += double(rand()) / RAND_MAX;
            sumA += double(rand()) / RAND_MAX;
            sumR += double(rand()) / RAND_MAX;
        }
        newW[l] = min(max((sumW - 6) * t + curW[l], 0.001), maxU);
        newA[l] = min(max((sumA - 6) * t + curA[l], minX), maxX);
        newR[l] = min(max((sumR - 6) * t + curR[l], 0.001), maxX);
    }
}

double TInitU::h(double deltaE, double t) {
    return exp(-deltaE / t);
}

double TInitU::E() {
    double res = 0;
    for (int i = 0; i < N; i++) {
        double d = U(newW, newA, newR, x[i]) - aU[i];
        res += d * d;
    }
    return res;
}

double TInitU::U(double *w, double *a, double *r, double x) {
    double res = 0;
    for (int j = 0; j < M; j++) {
        res += gaussian(w[j], a[j], r[j], x);
    }
    return res;
}

double TInitU::gaussian(double w, double a, double r, double x) {
	double d = x - a;
	return w * exp(-(d * d) / (r * r));
}
