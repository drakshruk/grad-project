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
#include <QMessageBox>
#include "dialog.h"
#include <QDebug>
#include "graphwidget.h"
#include "appdatamodel.h"
#include "imagecalculator.h"
#include "imageshowcasewidget.h"
#include "imageprocessor.h"
#include "loggerwidget.h"

namespace Ui {
class MainWindow;
}

class LoggerWidget;  // Forward declaration

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // Parameters
    int iRad = 3;
    double dSigma = 1.0;
    double dRedMax = 255.0;
    double dBlueMax = 255.0;
    bool edgeSelectionMode = false;

    // Matrices
    Matrix2D<double> matForTest;
    Matrix2D<double> imMat;
    Matrix2D<int> attMat;

    // Profile building points
    QVector<QPoint> profilePoints;
    void demonstrateBlurEffect();
    void demonstrateRadialBlurEffect();
    void demonstrateBlurEffectForRadiographicObject();
    void generate_wth();

    // Average time benchmarking
    void benchmarkConvolution();
    void benchmarkRefinement();

public slots:
    void on_actiongaussian_blur_triggered();
    void on_actiongaussian_edge_detection_triggered();
    void on_actionpcr_edge_detection_triggered();
    void on_actionopen_file_triggered();
    void on_actionsave_file_triggered();
    void on_showcaseWidget_clicked(const QPoint& imagePosition);
    void catch_ImageCalculator(const QImage& image1, const QImage& image2,
                               QString operation, bool newWindow, bool floatResult);
    void on_actiongradient_X_and_Y_triggered();
    void on_actiontest_triggered();
    void on_actiondraw_profile_triggered();
    void on_actionTwo_hollows_triggered();
    void on_actionTwo_hollows_big_triggered();
    void on_actionLaplacian_edge_detection_triggered();
    void on_imageUpdated(const QImage &newImage);
    void on_actionimage_calculator_triggered();
    void on_actionsharpen_triggered();
    void on_actiontest_002_triggered();

    // Slots for profile building
    void on_actionProfileBetweenPoints_triggered();
    void on_actionClearProfilePoints_triggered();
    void on_actionShowStatistics_triggered();
    void on_actionExportStatistics_triggered();

    // Slots for log window
    void on_actionShowLogWindow_triggered();
    void closeEvent(QCloseEvent *event);

private slots:
    void onImageClickedForProfile(const QPoint& imagePosition);

private:
    QImage drawProfileLineOnImage(const QImage& image, const QPoint& p1, const QPoint& p2, const QVector<double>& profile);
    // Helper methods
    void updateImageDisplay();
    void showStatisticsDialog(const ImageProcessor::ImageStatistics& stats);
    void showProfileDialog(const QVector<double>& profile, const QString& title);
    void resetEdgeSelection();

    // Member variables
    Ui::MainWindow  *ui;
    AppDataModel    *m_dataModel;
    QImage          curImage;
    QImage          im1;
    Dialog          *getDataDialog;
    GraphWidget     *GrWid;
    ImageCalculator *imCalculator;
    ImageShowcaseWidget *imageWidget;
    ImageShowcaseWidget *profileImageWidget;  // Separate widget for profile visualization
    LoggerWidget* m_loggerWidget = nullptr;

    // Edge selection (moved from global scope / peremeshcheno iz global'noy oblasti)
    QPoint mPos;
    QVector<QPoint> selectedEdge;
    QVector<QPoint> trueEdge;

    // Profile building state
    bool m_profileBuildingMode = false;
    QPoint m_firstProfilePoint;
    QPoint m_secondProfilePoint;
};

#endif // MAINWINDOW_H