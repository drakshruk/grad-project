#-------------------------------------------------
#
# Project created by QtCreator 2024-06-27T08:45:50
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Grad
TEMPLATE = app


SOURCES += src/main.cpp\
        src/appdatamodel.cpp \
        src/dialog.cpp \
        src/graphwidget.cpp \
        src/imagecalculator.cpp \
        src/imageprocessor.cpp \
        src/imageshowcasewidget.cpp \
        src/mainwindow.cpp \
        src/qcustomplot.cpp \
        src/utility.cpp

HEADERS  += include/mainwindow.h \
    include/appdatamodel.h \
    include/dialog.h \
    include/graphwidget.h \
    include/imagecalculator.h \
    include/imageprocessingtypes.h \
    include/imageprocessor.h \
    include/imageshowcasewidget.h \
    include/matrixlambdas.h \
    include/qcustomplot.h

FORMS    += ui/mainwindow.ui \
    ui/dialog.ui \
    ui/graphwidget.ui \
    ui/imagecalculator.ui \
    ui/imageshowcasewidget.ui
