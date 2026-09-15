#-------------------------------------------------
#
# Project created by QtCreator 2024-06-27T08:45:50
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Grad
TEMPLATE = app

# MinGW-specific OpenMP configuration
win32-g++ {
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp
    LIBS += -fopenmp
}

# MSVC configuration (if needed)
win32-msvc* {
    QMAKE_CXXFLAGS += /openmp
}

# Linux configuration
unix:!macx {
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp
}

CONFIG += c++17

SOURCES += main.cpp\
        appdatamodel.cpp \
        dialog.cpp \
        graphwidget.cpp \
        imagecalculator.cpp \
        imageprocessor.cpp \
        imageshowcasewidget.cpp \
        loggerwidget.cpp \
        mainwindow.cpp \
        qcustomplot.cpp

HEADERS  += mainwindow.h \
    appdatamodel.h \
    dialog.h \
    graphwidget.h \
    imagecalculator.h \
    imageprocessingtypes.h \
    imageprocessor.h \
    imageshowcasewidget.h \
    loggerwidget.h \
    matrixlambdas.h \
    qcustomplot.h

FORMS    += mainwindow.ui \
    dialog.ui \
    graphwidget.ui \
    imagecalculator.ui \
    imageshowcasewidget.ui \
    loggerwidget.ui

# Optimization flags
QMAKE_CXXFLAGS_RELEASE += -O3 -march=native
QMAKE_CXXFLAGS_DEBUG += -O0 -g