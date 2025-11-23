#-------------------------------------------------
#
# Project created by QtCreator 2024-06-27T08:45:50
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Grad
TEMPLATE = app


SOURCES += main.cpp\
        appdatamodel.cpp \
        dialog.cpp \
        graphwidget.cpp \
        imagecalculator.cpp \
        imageprocessor.cpp \
        imageshowcasewidget.cpp \
        mainwindow.cpp \
        qcustomplot.cpp \
        utility.cpp

HEADERS  += mainwindow.h \
    appdatamodel.h \
    dialog.h \
    graphwidget.h \
    imagecalculator.h \
    imageprocessingtypes.h \
    imageprocessor.h \
    imageshowcasewidget.h \
    matrixlambdas.h \
    qcustomplot.h

FORMS    += mainwindow.ui \
    dialog.ui \
    graphwidget.ui \
    imagecalculator.ui \
    imageshowcasewidget.ui
