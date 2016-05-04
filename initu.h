#ifndef INITU_H
#define INITU_H
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>

using namespace std;
class TInitU {
public:
    TInitU(int N, int DIM, double *x, double *aU, int maxIter, double eps);
    ~TInitU();
    double get(double* x);
private:
    void initRange();
    double findMinE();
    double findMinEm();
    void addGaussian();
    void firstGaussian();
    double E();
    void g();
    double h(double deltaE, double t);
    double T(double c, double T0);
    double gaussian(double w, double *a, double r, double *x);
    double U(double *w,double *a, double *r, double *x);
    int N, M, DIM;
    double deltaE;
    int maxIter;
    double *x;
    double *aU;

    double *newW;
    double *newA;
    double *newR;
    double *curW;
    double *curA;
    double *curR;

    double *minW;
    double *minA;
    double *minR;

    double eps;
    double t;
    double fromW, toW, fromR, toR;
    double *fromA, *toA;
    double minE;

    random_device rd;
    mt19937 gen;
};
#endif
