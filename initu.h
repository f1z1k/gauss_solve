#ifndef INITU_H
#define UNITU_H
 
class TInitU {
public:
    TInitU(int N, double *x, double *aU, int maxIter, double eps);
    ~TInitU();
    double get(double x);
private:
    double E();
    double findMinE();
    double findMinEm();
    void addGaussian();
    void firstGaussian();
    void g();
    double h(double deltaE, double t);
    double T(double c, double T0);
    double gaussian(double w, double a, double r, double x);
    double U(double *w,double *a, double *r, double x);
    int N, M;
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
    double minU, maxU, minX, maxX;
    double minE;
};

#endif
