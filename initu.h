#ifndef INITU_H
#define UNITU_H
 
class TInitU {
public:
    TInitU(int N, double *x, double *aU, int maxIter);
    ~TInitU();
    double get(double x);
private:
    double E();
    double findMinE();
    double findMinEm();
    void addWAR();
    void newWAR();
    double T(double c, double T0);
    double gaussian(double w, double a, double r, double x);
    double U(double *w,double *a, double *r, double x);
    int N, M;
    double deltaE;
    int maxIter;
    double *x;
    double *aU;

    double *prevW;
    double *prevA;
    double *prevR;
    double *curW;
    double *curA;
    double *curR;

    double *gradW;
    double *gradA;
    double *gradR;
};

#endif
