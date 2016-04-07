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

private:
    Ui::Widget *ui;
    int N, M;
    void plotAnalyze(double a, double b, double plotDelta);
    void plotResult(int N, int M, double a, double b, double plotDelta);
};

#endif // WIDGET_H
