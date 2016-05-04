#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include "initu.h"

namespace Ui {
class Widget;
}

class Widget : public QWidget
{
    Q_OBJECT

public:
    explicit Widget(QWidget *parent = 0);
    ~Widget();

private slots:
    void on_spinBox_valueChanged(int arg1);

    void on_spinBox_2_valueChanged(int arg1);

    void on_doubleSpinBox_valueChanged(double arg1);

    void on_pushButton_clicked();

    void on_doubleSpinBox_2_valueChanged(double arg1);

private:
    Ui::Widget *ui;
    int N, M;
    int DIM;
    double a, b;
    double plotDelta;
    double eps;
    double t;
    void plotMain();
    void plotAnalyze();
    void plotResult();
};

#endif // WIDGET_H
