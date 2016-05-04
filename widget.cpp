#include "widget.h"
#include "ui_widget.h"
#include "QDebug"

double U(double* x) {
    return sin(M_PI * x[0]) * (1 - exp(-M_PI * M_PI * x[1])) / (M_PI * M_PI);
}

double f(double* x) {
    return sin(M_PI * x[0]);
}

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);
    ui->plot->xAxis->setLabel("x");
    ui->plot->yAxis->setLabel("U");
    // set axes ranges, so we see all data:
    DIM = 2;
    N = 3;
    M = 1;
    a = 0;
    b = 1;
    t = 0;
    eps = 1e-10;
    plotDelta = .001;
    ui->plot->addGraph();
    ui->plot->addGraph();
    plotMain();
}

Widget::~Widget()
{
    delete ui;
}

void Widget::plotMain() {
    ui->plot->xAxis->setRange(a, b);
    ui->plot->yAxis->setRange(-1, 1);
    plotAnalyze();
    plotResult();
    ui->plot->replot();
}

void Widget::plotAnalyze() {
    int plotN = int((b - a) / plotDelta) + 1;
    QVector<double> plotX(plotN);
    QVector<double> plotU(plotN);
    double v[2];
    v[1] = t;
    int i = 0;
    for (double x = a; x <= b; x += plotDelta) {
        v[0] = x;
        plotX[i] = x;
        plotU[i] = U(v);
        i++;
    }

    ui->plot->graph(0)->setData(plotX, plotU);
    ui->plot->graph(0)->setPen(QPen(Qt::red));
}


void Widget::plotResult() {
    int plotN = int((b - a) / plotDelta) + 1;
    QVector<double> plotX(plotN);
    QVector<double> plotU(plotN);
    int i = 0;
    double x[1000];
    double testU[1000];
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> randX(a, b);
    uniform_real_distribution<> randT(0, t);

    for (int i = 0; i < N / 4; i++) {
        x[DIM * i] = randX(gen);
        x[DIM * (i + N / 4)] = a;
        x[DIM * (i + 2 * N / 4)] = b;

        x[DIM * i + 1] = 0;
        x[DIM * (i + N / 4) + 1] = randT(gen);
        x[DIM * (i + 2 * N / 4) + 1] = randT(gen);

        testU[i] = 0;
        testU[i + N / 4] = 0;
        testU[i + 2 * N / 4] = 0;
    }

    for (int i = 3 * N / 4; i < N; i++) {
        x[DIM * i] = randX(gen);
        x[DIM * i + 1] = randT(gen);
        testU[i] = f(x + DIM * i);
    }

    TInitU u(N, 2, x, testU, M, 1e-9);

    double v[2];
    v[1] = t;
    for (double x = a; x <= b; x += plotDelta) {
        plotX[i] = x;
        v[0] = x;
        double tmp = u.get(v);
        plotU[i] = tmp;
        i++;
    }

    ui->plot->graph(1)->setData(plotX, plotU);
    ui->plot->graph(1)->setPen(QPen(Qt::blue));
    return;
}



void Widget::on_spinBox_valueChanged(int arg1)
{
    N = arg1;
}

void Widget::on_spinBox_2_valueChanged(int arg1)
{
    M = arg1;
}

void Widget::on_doubleSpinBox_valueChanged(double arg1)
{
   eps = arg1;
}

void Widget::on_pushButton_clicked()
{
    plotMain();
}

void Widget::on_doubleSpinBox_2_valueChanged(double arg1)
{
    t = arg1;
}
