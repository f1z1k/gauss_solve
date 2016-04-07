#include "widget.h"
#include "ui_widget.h"
#include "QDebug"

double U(double x) {
    return 100 * sin(x);
}

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);
    double a = 0, b = 3;
    ui->plot->xAxis->setLabel("x");
    ui->plot->yAxis->setLabel("U");
    // set axes ranges, so we see all data:
    ui->plot->xAxis->setRange(a, b);
    ui->plot->yAxis->setRange(-1, b);
    N = 1000;
    M = 10;
    plotAnalyze(a, b, 0.001);
    plotResult(N, M, a, b, .001);
}

Widget::~Widget()
{
    delete ui;
}

void Widget::plotAnalyze(double a, double b, double plotDelta) {
    int plotN = int((b - a) / plotDelta) + 1;
    QVector<double> plotX(plotN);
    QVector<double> plotU(plotN);
    int i = 0;
    for (double x = a; x <= b; x += plotDelta) {
        plotX[i] = x;
        plotU[i] = U(x);
        i++;
    }

    ui->plot->addGraph();
    ui->plot->graph(0)->setData(plotX, plotU);
}

void Widget::plotResult(int N, int M, double a, double b, double plotDelta) {
    int plotN = int((b - a) / plotDelta) + 1;
    QVector<double> plotX(plotN);
    QVector<double> plotU(plotN);
    int i = 0;
    double x[10000];
    double testU[10000];
    double delta = (b - a) / N;
    for (int i = 0; i < N; i++) {
        x[i] = i * delta;
        testU[i] = U(x[i]);
    }

    TInitU u(N, x, testU, M);

    for (double x = a; x <= b; x += plotDelta) {
        plotX[i] = x;
        double tmp = u.get(x);
        plotU[i] = tmp;
        i++;
    }

    ui->plot->addGraph();
    ui->plot->graph(1)->setData(plotX, plotU);
}



void Widget::on_spinBox_valueChanged(int arg1)
{
    N = arg1;
    plotResult(N, M, 0, 1, 0.001);
}

void Widget::on_spinBox_2_valueChanged(int arg1)
{
   M = arg1;
   plotResult(N, M, 0, 1, 0.001);

}
