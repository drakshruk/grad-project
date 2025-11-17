#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QImage>
#include <QLabel>
#include <QPixmap>
#include <QList>
#include <cmath>
#include <iostream>
#include <QFile>
#include <QTextStream>
#include <QMouseEvent>
#include <QFileDialog>
#include "dialog.h"
#include <QDebug>
#include "graphwidget.h"
#include "appdatamodel.h"
#include "imagecalculator.h"
#include "imageshowcasewidget.h"
#include "imageprocessor.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    int iRad;
    double dSigma, dRedMax, dBlueMax;
    bool edgeSelectionMode = false;
    Matrix2D<double> matForTest;
    Matrix2D<double> imMat;
    Matrix2D<int> attMat;

public slots:

    void on_actiongaussian_blur_triggered();

    void on_actiongaussian_edge_detection_triggered();

    void on_actionpcr_edge_detection_triggered();

    void on_actionopen_file_triggered();

    void on_actionsave_file_triggered();

    void mousePressEvent(QMouseEvent *e);

    void catch_ImageCalculator(const QImage& image1, const QImage& image2, QString operation, bool newWindow, bool floatResult);

    void on_actiongradient_X_and_Y_triggered();

    void on_actiontest_triggered();

    void on_actiondraw_profile_triggered();

    void on_actionTwo_hollows_triggered();

    void on_actionTwo_hollows_big_triggered();

    void on_actionLaplacian_edge_detection_triggered();

    void on_imageUpdated(const QImage &newImage);

    void on_actionimage_calculator_triggered();

    void on_actionsharpen_triggered();

private:
    Ui::MainWindow  *ui;
    AppDataModel    *m_dataModel;
    QPixmap         *pixmap;
    QImage          curImage, im1;
    Dialog          *getDataDialog;
    GraphWidget     *GrWid;
    ImageCalculator *imCalculator;
    ImageShowcaseWidget *imageWidget;
};

#endif // MAINWINDOW_H
