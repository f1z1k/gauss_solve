#include <iostream>

#include <cmath>
#include <cstring>
#include <QDebug>

#include "initu.h"

using namespace std;

// public
TInitU::TInitU(int N, int DIM, double *x, double *aU, int maxIter, double eps)
    :N(N), DIM(DIM), maxIter(maxIter), M(maxIter), eps(eps), gen(rd())
{
    this->x     = new double [N * DIM];
    this->aU    = new double [N * DIM];
    this->newW  = new double [maxIter];
    this->newA  = new double [maxIter * DIM];
    this->newR  = new double [maxIter];
    this->curW  = new double [maxIter];
    this->curA  = new double [maxIter * DIM];
    this->curR  = new double [maxIter];
    this->minW  = new double [maxIter];
    this->minA  = new double [maxIter * DIM];
    this->minR  = new double [maxIter];
    this->fromA = new double [maxIter * DIM];
    this->toA   = new double [maxIter * DIM];
    memcpy(this->x, x, sizeof(double) * N * DIM);
    memcpy(this->aU, aU, sizeof(double) * N * DIM);

    initRange();
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
    delete[] this->fromA;
    delete[] this->toA;
}

double TInitU::get(double* x) {
    return U(minW, minA, minR, x);
}

// private


void TInitU::initRange() {

    fromW = toW = aU[0];
    fromR = toR = eps;
    memcpy(fromA, x, sizeof(double) * DIM);
    memcpy(toA, x, sizeof(double) * DIM);
    for (int i = 1; i < N; i++) {
        fromW = min(fromW, aU[i]);
        toW = max(toW, aU[i]);
        for (int l = 0; l < DIM; l++) {
            fromA[l] = min(fromA[l], x[i * DIM + l]);
            toA[l]   = max(toA[l], x[i * DIM + l]);

            toR = max(toR, abs(x[i * DIM + l]));
        }
    }
}

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
        qDebug() << minW[l] << " " << minR[l];
        for (int j = 0; j < DIM; j++) {
            qDebug() << minA[l * DIM + j];
        }
    }
    return minE;
}

double TInitU::findMinEm() {
    double curE = E();
    t = T(1, 1000 * toR);
    while (t > eps) {
        t = T(0.99, t);
        if (curE < minE) {
            minE = curE;
            memcpy(minW, curW, sizeof(double) * M);
            memcpy(minA, curA, sizeof(double) * M * DIM);
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
    uniform_real_distribution<> w(fromW, toW);
    uniform_real_distribution<> r(fromR, toR);
    minW[0] = curW[0] = newW[0] = w(gen);
    minR[0] = curR[0] = newR[0] = r(gen);
    for (int l = 0; l < DIM; l++) {
        uniform_real_distribution<> a(fromA[l], toA[l]);
        minA[l] = curA[l] = newA[l] = a(gen);
    }
}

void TInitU::addGaussian() {
    double k = 2;
    minW[M - 1] = curW[M - 1] = newW[M - 1] = eps;
    minR[M - 1] = curR[M - 1] = newR[M - 1] = k * fromR;
    for (int l = 0; l < DIM; l++) {
        int i = (M - 1) * DIM + l;
        minA[i] = curA[i] = newA[i] = 0.5 * k * (fromA[l] + toA[l]);
    }
}

double TInitU::T(double c, double prevT) {
    return c * prevT;
}

double normalRand(double m, double r) {
    static int status = 0;
    static double r1;
    static double r2;
    static double n1;
    static double n2;
    status++;
    if ((status & 1)) {
        r1 = double(rand()) / RAND_MAX;
        r2 = double(rand()) / RAND_MAX;
        n1 = -2 * log(r1) * sin(2 * M_PI * r2);
        n2 = -2 * log(r1) * cos(2 * M_PI * r2);
        return m + r * n1;
    }
    return m + r * n2;
}

void TInitU::g() {
    int k = 1;

    for (int l = 0; l < M; l++) {
        normal_distribution<> normalW(curW[l], t);
        normal_distribution<> normalR(curR[l], t);
        newW[l] = min(max(normalW(gen), k * fromW), k * toW);
        newR[l] = min(max(normalR(gen), fromR), k * toR);
        for (int j = 0; j < DIM; j++) {
            normal_distribution<> normalA(curA[l * DIM + j], t);
            newA[l * DIM + j] = min(max(normalA(gen),  k * fromA[j]), k * toA[j]);
        }
    }
}

double TInitU::h(double deltaE, double t) {
    return exp(-deltaE / t);
}

double TInitU::E() {
    double res = 0;
    for (int i = 0; i < N / 4; i++) {
        double d = U(newW, newA, newR, x + i * DIM) - aU[i];
        res += d * d;
    }
    for (int i = N / 4; i < 2 * N / 4; i++) {
        double d = U(newW, newA, newR, x + i * DIM) - aU[i];
        res += d * d;
    }
    for (int i = 2 * N / 4; i < 3 * N / 4; i++) {
        double d = U(newW, newA, newR, x + i * DIM) - aU[i];
        res += d * d;
    }
    for (int i = 3 * N / 4; i < N; i++) {
        double* v = x + i * DIM;
        double d = 0;
        for (int j = 0; j < M; j++) {
            double w = newW[j];
            double *a = newA + j * DIM;
            double r = newR[j];
            double k1 = -2 * (v[1] - a[1]) / (r * r);
            double k2 = (-2 * r * r + 4 * a[0] * a[0] - 8 * a[0] * v[0] + 4 * v[0] * v[0]) / (r * r * r * r);
            d += gaussian((k1 - k2) * w, a, r, v);
        }
        d -= aU[i];
        res += d * d;
    }
    return res;
}

double TInitU::U(double *w, double *a, double *r, double *x) {
    double res = 0;
    for (int j = 0; j < M; j++) {
        res += gaussian(w[j], a + j * DIM, r[j], x);
    }
    return res;
}

double TInitU::gaussian(double w, double* a, double r, double* x) {
    double sum = 0;
    for (int j = 0; j < DIM; j++) {
        double d = x[j] - a[j];
        sum += d * d;
    }
    return w * exp(-sum / (r * r));
}
