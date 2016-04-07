#-------------------------------------------------
#
# Project created by QtCreator 2016-03-28T23:42:19
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = solve_gauss
TEMPLATE = app


SOURCES += main.cpp\
        widget.cpp \
    initu.cpp \
    qcustomplot.cpp

HEADERS  += widget.h \
    initu.h \
    qcustomplot.h

FORMS    += widget.ui
