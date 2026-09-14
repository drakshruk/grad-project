#include "mainwindow.h"
#include "ui_mainwindow.h"

QPoint mPos;
QVector<QPoint> selectedEdge;
QVector<QPoint> trueEdge;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_dataModel(new AppDataModel(this))
{
    ui->setupUi(this);

    // Initialize widgets
    imCalculator = new ImageCalculator();
    imageWidget = new ImageShowcaseWidget();
    profileImageWidget = new ImageShowcaseWidget();
    getDataDialog = new Dialog();
    GrWid = new GraphWidget();

    // Set up connections
    connect(m_dataModel, &AppDataModel::currentImageChanged,
            this, &MainWindow::on_imageUpdated);
    connect(imCalculator, &ImageCalculator::throw_imageCalculator,
            this, &MainWindow::catch_ImageCalculator);
    connect(imageWidget, &ImageShowcaseWidget::imageClicked,
            this, &MainWindow::on_showcaseWidget_clicked);

    // Initialize parameters
    m_dataModel->setBlueMax(255);
    m_dataModel->setRedMax(255);
    m_dataModel->setRadius(10);
    m_dataModel->setSigma(5);

    // Set window properties
    this->setFixedSize(550, 50);

    // Update parameters from model
    iRad = m_dataModel->radius();
    dSigma = m_dataModel->sigma();
    dRedMax = m_dataModel->redMax();
    dBlueMax = m_dataModel->blueMax();
}

MainWindow::~MainWindow()
{
    delete ui;
    delete profileImageWidget;
}

// ============================================================================
// File Operations
// ============================================================================

void MainWindow::on_actionopen_file_triggered()
{
    QFileDialog* openImDialog = new QFileDialog(this);
    openImDialog->setFileMode(QFileDialog::AnyFile);
    openImDialog->setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));

    QString fileName = openImDialog->getOpenFileName();
    if (fileName.isEmpty()) return;

    QImage newIm(fileName);
    if (newIm.isNull()) {
        QMessageBox::warning(this, "Error", "Failed to load image: " + fileName);
        return;
    }

    m_dataModel->setCurrentImage(newIm);
    m_dataModel->setOriginalImage(newIm);
    edgeSelectionMode = false;
    m_profileBuildingMode = false;
}

void MainWindow::on_actionsave_file_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image to save");
        return;
    }

    QFileDialog* saveImDialog = new QFileDialog(this);
    saveImDialog->setFileMode(QFileDialog::AnyFile);
    saveImDialog->setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));
    QString saveFileName = saveImDialog->getSaveFileName();

    if (!saveFileName.isEmpty()) {
        curImage.save(saveFileName);
    }
}

// ============================================================================
// Image Processing Operations
// ============================================================================

void MainWindow::on_actiongaussian_blur_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    getDataDialog->setWindowTitle("Write sigma coefficient value");
    getDataDialog->setPlaceholderText("Sigma value");
    getDataDialog->setDefaultValue(5.0);

    if (getDataDialog->exec()) {
        double sigma = getDataDialog->getValue();
        m_dataModel->setSigma(sigma);
        m_dataModel->setRadius(8 * m_dataModel->sigma());

        curImage = ImageProcessor::convImage(curImage,
            ImageProcessor::getGauss(m_dataModel->radius(), m_dataModel->radius(), m_dataModel->sigma()));

        updateImageDisplay();
        im1 = curImage;

        qDebug() << "Gaussian blur applied with sigma =" << sigma;
    }
}

void MainWindow::on_actiongaussian_edge_detection_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    getDataDialog->setWindowTitle("Write sigma coefficient value");
    getDataDialog->setPlaceholderText("Sigma value");
    getDataDialog->setDefaultValue(5.0);

    if (getDataDialog->exec()) {
        m_dataModel->setSigma(getDataDialog->getValue());
        m_dataModel->setRadius(8 * m_dataModel->sigma());

        curImage = ImageProcessor::gaussianEdgeDetection(
            ImageProcessor::fromGrayImage(curImage),
            m_dataModel->sigma(),
            m_dataModel->radius(),
            attMat);

        updateImageDisplay();
        edgeSelectionMode = true;

        qDebug() << "Gaussian edge detection applied";
    }
}

void MainWindow::on_actionLaplacian_edge_detection_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    getDataDialog->setWindowTitle("Write sigma coefficient value");
    getDataDialog->setPlaceholderText("Sigma value");
    getDataDialog->setDefaultValue(3.0);

    if (getDataDialog->exec()) {
        m_dataModel->setSigma(getDataDialog->getValue());
        m_dataModel->setRadius(3 * m_dataModel->sigma());

        Matrix2D<double> mat1 = ImageProcessor::fromGrayImage(curImage);
        Matrix2D<double> mat2 = ImageProcessor::convMat(mat1,
            ImageProcessor::getLapl(m_dataModel->radius(), m_dataModel->radius(), m_dataModel->sigma()));

        matForTest = mat2;
        mat1 = ImageProcessor::elementWiseOperation(mat1, mat2, MatrixLambdas::Subtract<double>{});
        matForTest = mat1;
        mat1 = ImageProcessor::findEdges(mat1, 3);
        curImage = ImageProcessor::toGrayImage(mat1);

        updateImageDisplay();
        edgeSelectionMode = true;

        qDebug() << "Laplacian edge detection applied";
    }
}

void MainWindow::on_actionsharpen_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    Matrix2D<double> kernel = {{-1, -1, -1},
                               {-1,  9, -1},
                               {-1, -1, -1}};
    curImage = ImageProcessor::convImage(curImage, kernel);

    updateImageDisplay();
    im1 = curImage;

    qDebug() << "Sharpen filter applied";
}

void MainWindow::on_actiongradient_X_and_Y_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    Matrix2D<double> xGrad = ImageProcessor::fromGrayImage(curImage);
    xGrad = ImageProcessor::convMat(xGrad, ImageProcessor::getXGradCore(iRad, iRad, dSigma));
    ImageProcessor::toBlueRedImage(xGrad, 255., 255.).save("XGradImage.png");

    Matrix2D<double> yGrad = ImageProcessor::fromGrayImage(curImage);
    yGrad = ImageProcessor::convMat(yGrad, ImageProcessor::getYGradCore(iRad, iRad, dSigma));
    ImageProcessor::toBlueRedImage(yGrad, 255., 255.).save("YGradImage.png");

    qDebug() << "Gradient images saved to XGradImage.png and YGradImage.png";
}

// ============================================================================
// Sample Images
// ============================================================================

void MainWindow::on_actionTwo_hollows_triggered()
{
    QImage setImage = ImageProcessor::sampleTwoHollows();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);
    edgeSelectionMode = false;
    m_profileBuildingMode = false;


    // First circle (bottom) / Pervaya okruzhnost' (snizu)
    int centerX1 = 100, centerY1 = 135, radius1 = 40;
    // Second circle (top) / Vtoraya okruzhnost' (sverkhu)
    int centerX2 = 100, centerY2 = 65, radius2 = 30;

    for(int i = 0; i < 200; i++){
        for(int j = 0; j < 200; j++){
            // Check first circle / Proveryaem pervuyu okruzhnost'
            int dx1 = i - centerX1;
            int dy1 = j - centerY1;
            int distSq1 = dx1*dx1 + dy1*dy1;

            if(abs(distSq1 - radius1 * radius1) <= 1){
                trueEdge.push_back({i,j});
            }

            // Check second circle / Proveryaem vtoruyu okruzhnost'
            int dx2 = i - centerX2;
            int dy2 = j - centerY2;
            int distSq2 = dx2*dx2 + dy2*dy2;

            if(abs(distSq2 - radius2 * radius2) <= 1){
                trueEdge.push_back({i,j});
            }
        }
    }

}

void MainWindow::on_actionTwo_hollows_big_triggered()
{
    QImage setImage = ImageProcessor::sampleTwoHollowsBig();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);
    edgeSelectionMode = false;
    m_profileBuildingMode = false;
}

// ============================================================================
// Profile Operations
// ============================================================================

void MainWindow::on_actiondraw_profile_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    QVector<double> x, y;
    int centerX = curImage.width() / 2;

    for (int i = 0; i < curImage.height(); i++) {
        x.push_back(i);
        y.push_back(qGray(curImage.pixel(centerX, i)));
    }

    GrWid->plotGraph(x, y);
    GrWid->setAxisLabels("Y coordinate", "Intensity");
    GrWid->setTitle(QString("Horizontal Profile at x = %1").arg(centerX));
    GrWid->show();
}

void MainWindow::on_actionProfileBetweenPoints_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    if (profilePoints.size() < 2) {
        m_profileBuildingMode = true;
        profilePoints.clear();/*
        QMessageBox::information(this, "Profile Building",
            "Click two points on the image to build a profile.\n"
            "The profile will show intensity values along the line between the points.");*/
    } else {
        // Build profile between the two selected points
        QVector<double> profile = ImageProcessor::buildProfileBetweenPoints(
            profilePoints[0], profilePoints[1], curImage);

        if (profile.isEmpty()) {
            QMessageBox::warning(this, "Error", "Failed to build profile");
            return;
        }

        // Display profile
        QVector<double> x(profile.size());
        for (int i = 0; i < profile.size(); i++) {
            x[i] = i;
        }

        GrWid->plotGraph(x, profile);
        GrWid->setAxisLabels("Distance along line (pixels)", "Intensity");
        GrWid->setTitle(QString("Profile from (%1,%2) to (%3,%4)")
            .arg(profilePoints[0].x()).arg(profilePoints[0].y())
            .arg(profilePoints[1].x()).arg(profilePoints[1].y()));
        GrWid->show();

        // Also show statistical summary
        double minVal = *std::min_element(profile.begin(), profile.end());
        double maxVal = *std::max_element(profile.begin(), profile.end());
        double sum = 0.0;
        for (double v : profile) sum += v;
        double mean = sum / profile.size();

        qDebug() << "Profile Statistics:";
        qDebug() << "  Length:" << profile.size() << "pixels";
        qDebug() << "  Min:" << minVal;
        qDebug() << "  Max:" << maxVal;
        qDebug() << "  Mean:" << mean;

        // Reset for next profile
        profilePoints.clear();
        m_profileBuildingMode = false;
    }
}

void MainWindow::on_actionClearProfilePoints_triggered()
{
    profilePoints.clear();
    m_profileBuildingMode = false;
//    QMessageBox::information(this, "Profile Building", "Profile points cleared.");
}

void MainWindow::onImageClickedForProfile(const QPoint& imagePosition)
{
    if (!m_profileBuildingMode) return;

    profilePoints.append(imagePosition);

    if (profilePoints.size() == 1) {/*
        QMessageBox::information(this, "Profile Building",
            QString("First point selected at (%1,%2). Click second point.")
            .arg(imagePosition.x()).arg(imagePosition.y()));*/
    } else if (profilePoints.size() == 2) {
        // Build and display profile
        QVector<double> profile = ImageProcessor::buildProfileBetweenPoints(
            profilePoints[0], profilePoints[1], curImage);

        if (!profile.isEmpty()) {
            QVector<double> x(profile.size());
            for (int i = 0; i < profile.size(); i++) x[i] = i;

            GrWid->plotGraph(x, profile);
            GrWid->setAxisLabels("Distance along line (pixels)", "Intensity");
            GrWid->setTitle(QString("Profile from (%1,%2) to (%3,%4)")
                .arg(profilePoints[0].x()).arg(profilePoints[0].y())
                .arg(profilePoints[1].x()).arg(profilePoints[1].y()));
            GrWid->show();
        }

        profilePoints.clear();
        m_profileBuildingMode = false;
    }
}

// ============================================================================
// Statistics Operations
// ============================================================================

void MainWindow::on_actionShowStatistics_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    auto stats = ImageProcessor::computeImageStatistics(curImage, attMat);
    showStatisticsDialog(stats);
}

void MainWindow::on_actionExportStatistics_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    QString fileName = QFileDialog::getSaveFileName(this, "Save Statistics",
        "", "Text Files (*.txt);;All Files (*)");

    if (fileName.isEmpty()) return;

    auto stats = ImageProcessor::computeImageStatistics(curImage, attMat);
    if (ImageProcessor::saveStatisticsToFile(stats, fileName)) {
        QMessageBox::information(this, "Success", "Statistics saved to " + fileName);
    } else {
        QMessageBox::warning(this, "Error", "Failed to save statistics");
    }
}

void MainWindow::on_actionExportHistogram_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    QString fileName = QFileDialog::getSaveFileName(this, "Export Histogram",
        "", "CSV Files (*.csv);;All Files (*)");

    if (fileName.isEmpty()) return;

    auto stats = ImageProcessor::computeImageStatistics(curImage, attMat);
    if (ImageProcessor::exportHistogramToCSV(stats, fileName)) {
        QMessageBox::information(this, "Success", "Histogram exported to " + fileName);
    } else {
        QMessageBox::warning(this, "Error", "Failed to export histogram");
    }
}

// ============================================================================
// Edge Selection and Refinement
// ============================================================================

void MainWindow::on_showcaseWidget_clicked(const QPoint &imagePosition)
{
    mPos.setX(imagePosition.x());
    mPos.setY(imagePosition.y());

    // Handle profile building mode
    if (m_profileBuildingMode) {
        onImageClickedForProfile(imagePosition);
        return;
    }

    // Handle edge selection mode
    if (edgeSelectionMode && !attMat.empty()) {
        resetEdgeSelection();

        selectedEdge.clear();
        for(int i = 0; i < attMat.size(); i++) {
            for(int j = 0; j < attMat[0].size(); j++) {
                if(attMat[i][j] == static_cast<int>(attribute::isEdge))
                    selectedEdge.push_back({i,j});
            }
        }
//        ImageProcessor::selectEdge(mPos, attMat, selectedEdge);

        qDebug() << "Selected edge size =" << selectedEdge.size();

        m_dataModel->setSelectedEdge(selectedEdge);

        // Highlight selected edge in cyan
        for (int i = 0; i < selectedEdge.size(); i++) {
            curImage.setPixel(selectedEdge[i], qRgb(0, 255, 255));
        }

        updateImageDisplay();
    }
}

void MainWindow::resetEdgeSelection()
{
    if (attMat.empty()) return;

    for (int i = 0; i < curImage.width(); i++) {
        for (int j = 0; j < curImage.height(); j++) {
            if (attMat[i][j] == static_cast<int>(attribute::isSelectedEdge)) {
                attMat[i][j] = static_cast<int>(attribute::isEdge);
                curImage.setPixel(i, j, qRgb(255, 255, 255));
            }
        }
    }
}

void MainWindow::on_actiontest_triggered()
{
    if (!edgeSelectionMode) {
        QMessageBox::warning(this, "Error", "Edge selection mode is off. Run edge detection first.");
        return;
    }

    // Collect all selected edge points
    QVector<QPoint> selectedEdgePoints = selectedEdge;
    for (int i = 0; i < (int)attMat.size(); i++) {
        for (int j = 0; j < (int)attMat[i].size(); j++) {
            if (attMat[i][j] == static_cast<int>(attribute::isSelectedEdge)) {
                selectedEdgePoints.append(QPoint(i, j));
            }
        }
    }

    if (selectedEdgePoints.isEmpty()) {
        QMessageBox::warning(this, "Error", "No edge points selected");
        return;
    }

    qDebug() << "Refining" << selectedEdgePoints.size() << "edge points";

    // Create white background image for visualization
    QImage whiteImage = m_dataModel->originalImage();
//    for (int i = 0; i < whiteImage.width(); i++) {
//        for (int j = 0; j < whiteImage.height(); j++) {
//            if (qGray(whiteImage.pixel(i, j)) != 0) {
//                whiteImage.setPixel(i, j, qRgb(255, 255, 255));
//            } else {
//                whiteImage.setPixel(i, j, qRgb(0, 0, 0));
//            }
//        }
//    }

    // Image dimensions
    int NX = m_dataModel->originalImage().width();
    int NY = m_dataModel->originalImage().height();

    // Gaussian blur parameters
    double sigma0 = 8.0;
    double sigma01 = sqrt(sigma0 * sigma0 + 4.0 * 4.0);
    int NS = 2 * int(4 * sigma01 + 0.5);

    // Two different sigma values for edge detection
    double sigma1 = 4.0;
    double sigma2 = 8.0;
    double sigma_max = max(sigma1, sigma2);

    // Prepare matrices
    Matrix2D<double> A0 = ImageProcessor::fromGrayImage(m_dataModel->originalImage());
    Matrix2D<double> B1 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma0));
    Matrix2D<double> B01 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma01));
    Matrix2D<double> A = B1;

    // Update kernel size for max sigma
    NS = 2 * int(4 * sigma_max + 0.5);

    // Calculate gradients
    Matrix2D<double> B = A;
    Matrix2D<double> GxW1 = ImageProcessor::convMat(B, ImageProcessor::getXGradCore(NS, NS, sigma1));
    Matrix2D<double> GyW1 = ImageProcessor::convMat(B, ImageProcessor::getYGradCore(NS, NS, sigma1));

    // Profile parameters
    int n_sigma = 10;
    int n_myu = 200;
    double sigma_myu = min(sigma1, sigma2);
    int NN = 100;
    int otstup = 15;

    RefinementParameters params;
    params.A = A;
    params.B01 = B01;
    params.NN = NN;
    params.NX = NX;
    params.NY = NY;
    params.n_myu = n_myu;
    params.n_sigma = n_sigma;
    params.otstup = otstup;
    params.sigma0 = sigma0;
    params.sigma01 = sigma01;
    params.sigma1 = sigma1;
    params.sigma2 = sigma2;

    QVector<QPointF> refinedPositions;

    // Process each selected point
    for (int pointIdx = 0; pointIdx < selectedEdgePoints.size(); pointIdx++) {
        int n0 = selectedEdgePoints[pointIdx].x();
        int m0 = selectedEdgePoints[pointIdx].y();

        double ex = GxW1[n0][m0];
        double ey = GyW1[n0][m0];
        double gradMag = sqrt(ex * ex + ey * ey);

        if (gradMag > 1e-6) {
            ex /= gradMag;
            ey /= gradMag;
        } else {
            ex = 1.0;
            ey = 0.0;
        }

        params.ex = ex;
        params.ey = ey;

        RefinementResult result = ImageProcessor::refineSinglePoint(n0, m0, params);
        double n_new = result.refinedPosition.x();
        double m_new = result.refinedPosition.y();
        refinedPositions.append(QPointF(n_new, m_new));

        // Visualization
        if (n0 >= 0 && n0 < whiteImage.width() && m0 >= 0 && m0 < whiteImage.height()) {
            whiteImage.setPixel(n0, m0, qRgb(0, 0, 255));  // Blue - gradient
        }

        if (n_new >= 0 && n_new < whiteImage.width() && m_new >= 0 && m_new < whiteImage.height()) {
            whiteImage.setPixel(n_new, m_new, qRgb(255, 0, 0));  // Red - refined
        }

        qDebug() << "Point" << pointIdx + 1 << "/" << selectedEdgePoints.size()
                 << "refined: (" << n0 << "," << m0 << ") -> ("
                 << n_new << "," << m_new << ") shift =" << result.FFF1_1;
    }

    for(int i = 0; i < trueEdge.size(); i++) {
        whiteImage.setPixel(trueEdge[i].x(), trueEdge[i].y(), qRgb(0, 255, 0));  // Green - true
    }

    // Statistics summary
    qDebug() << "\n=== REFINEMENT SUMMARY ===";
    qDebug() << "Total points processed:" << refinedPositions.size();

    double avgShiftX = 0.0, avgShiftY = 0.0;
    for (int i = 0; i < selectedEdgePoints.size(); i++) {
        double shiftX = refinedPositions[i].x() - selectedEdgePoints[i].x();
        double shiftY = refinedPositions[i].y() - selectedEdgePoints[i].y();
        avgShiftX += shiftX;
        avgShiftY += shiftY;
        qDebug() << "Point" << i+1 << ": (" << selectedEdgePoints[i].x() << "," << selectedEdgePoints[i].y()
                 << ") -> (" << refinedPositions[i].x() << "," << refinedPositions[i].y()
                 << ") shift (" << shiftX << "," << shiftY << ")";
    }

    if (refinedPositions.size() > 0) {
        avgShiftX /= refinedPositions.size();
        avgShiftY /= refinedPositions.size();
        qDebug() << "Average shift: (" << avgShiftX << "," << avgShiftY << ")";
    }

    // Show final result
    ImageShowcaseWidget* resultWidget = new ImageShowcaseWidget();
    resultWidget->setImage(whiteImage);
    resultWidget->setWindowTitle(QString("Edge Refinement Result - %1 points").arg(refinedPositions.size()));
    resultWidget->show();

    qDebug() << "Edge refinement completed";
}

void MainWindow::on_actiontest_002_triggered()
{
    // Generate test image
    QImage setImage = ImageProcessor::sampleTwoHollows();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);

    // Image dimensions
    int NX = setImage.width();
    int NY = setImage.height();

    // Gaussian blur parameters
    double sigma0 = 8.0;
    double sigma01 = sqrt(sigma0 * sigma0 + 4.0 * 4.0);
    int NS = 2 * int(4 * sigma01 + 0.5);
    qDebug() << "NS =" << NS;

    // Two different sigma values for edge detection
    double sigma1 = 4.0;
    double sigma2 = 8.0;
    double sigma_max = max(sigma1, sigma2);

    // Convert images to matrices and apply Gaussian blur
    Matrix2D<double> A0 = ImageProcessor::fromGrayImage(m_dataModel->currentImage());
    Matrix2D<double> B1 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma0));
    Matrix2D<double> B01 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma01));
    Matrix2D<double> A = B1;

    // Show profiles before and after convolution
    GraphWidget* grwid1 = new GraphWidget();
    QVector<double> profileA, profileA0, profileX;
    for (int i = 0; i < (int)A.size(); i++) {
        profileA.push_back(A[100][i]);
        profileA0.push_back(A0[100][i]);
        profileX.push_back(i);
    }
    grwid1->plotTwoGraphs(profileX, profileA, profileX, profileA0);
    grwid1->setWindowTitle("Profiles before and after convolution");
    grwid1->show();

    // Update kernel size for max sigma
    NS = 2 * int(4 * sigma_max + 0.5);
    qDebug() << "NS =" << NS;

    // Perform Gaussian edge detection
    curImage = ImageProcessor::gaussianEdgeDetection(
        ImageProcessor::fromGrayImage(curImage), sigma_max, NS, attMat);
    edgeSelectionMode = true;
    updateImageDisplay();

    // Select test point for refinement
    int n0 = 78, m0 = 83;

    // Calculate gradients
    Matrix2D<double> B = A;
    Matrix2D<double> GxW1 = ImageProcessor::convMat(B, ImageProcessor::getXGradCore(NS, NS, sigma1));
    Matrix2D<double> GyW1 = ImageProcessor::convMat(B, ImageProcessor::getYGradCore(NS, NS, sigma1));

    double ex = GxW1[n0][m0];
    double ey = GyW1[n0][m0];
    double gradMag = sqrt(ex * ex + ey * ey);

    if (gradMag > 1e-6) {
        ex /= gradMag;
        ey /= gradMag;
    } else {
        ex = 1.0;
        ey = 0.0;
    }

    qDebug() << "x grad at (n0,m0):" << ex;
    qDebug() << "y grad at (n0,m0):" << ey;

    qDebug() << "180./pi*acos(ex) = " << 180./pi*acos(ex) << "\n";
    qDebug() << "180./pi*acos(ex) = " << 180./pi*acos(ey) << "\n";


    //*****************************************
    // PROFILE CONSTRUCTION / POSTROENIE PROFILEY
    //*****************************************

    QVector<double> prof1, prof2;
    int n_sigma = 10;           // Number of sigma steps / Kolichestvo shagov sigma
    int n_myu = 200;             // Number of mu points / Kolichestvo tochek mu
    int x0 = n0, y0 = m0;        // Starting point / Nachalnaya tochka
    double sigma_myu = min(sigma1, sigma2);  // Sigma for profile / Sigma dlya profilya
    QVector<double> mu;

    // Generate mu values / Generiruem znacheniya mu
    for(int s = 0; s < n_myu; s++) {
        double val = -n_sigma*sigma_myu + s*(2*n_sigma*sigma_myu)/n_myu;
        mu.push_back(val);
    }

    // Build profiles along gradient direction / Stroim profili vdol napravleniya gradienta
    for(int s = 0; s < n_myu; s++) {
        double prof1Val = 0.0;
        double prof2Val = 0.0;
        for(int n = 0; n < NX; n++) {
            for(int m = 0; m < NY; m++) {
                double mult = 0.0;

                // Calculate distance from point along gradient / Vychislyaem rasstoyanie ot tochki vdol gradienta
                double dx = (x0 - n + mu[s]*ex);
                double dy = (y0 - m + mu[s]*ey);
                double rad = dx*dx + dy*dy;

                // Laplacian of Gaussian kernel / Yadro Laplaciana Gaussa
                mult = (rad/sigma1/sigma1 - 2.0) * exp(-1./2.*rad/sigma1/sigma1);

                prof1Val += A[n][m]*mult;
                prof2Val += B01[n][m]*mult;
            }
        }
        prof1.push_back(prof1Val);
        prof2.push_back(prof2Val);
    }

    //*****************************************

    // Show initial profiles / Pokazyvaem iskhodnye profili
    GraphWidget* grwid2 = new GraphWidget();
    QVector<double> x1, x2, y1, y2;
    for(int i = 0; i < n_myu; i++) {
        x1.push_back(i);
        x2.push_back(i);
        y1.push_back(prof1[i]);
        y2.push_back(prof2[i]);
    }
    grwid2->plotTwoGraphs(x1, y1, x2, y2);
    grwid2->setWindowTitle(QString("Profiles before aligning"));
    grwid2->show();

    // Find zero crossings / Nahodim perehody cherez nol
    App_Stats y2_stats;
    y2_stats.gather_stats(y2);
    int nmumax = max(y2_stats.x_max, y2_stats.x_min), nmumin = min(y2_stats.x_max, y2_stats.x_min);
    int N0 = n_myu;

    qDebug() << "nmumin = " << nmumin << "nmumax = " << nmumax;
    int nL_Zero = N0/2, nR_Zero = N0/2;
    for(int n = nmumin; n < nmumax; n++) {
        if(prof1[n]*prof1[n+1] < 0) nL_Zero = n;
        if(prof2[n]*prof2[n+1] < 0) nR_Zero = n;
    }

    int N_Zero = (nL_Zero + nR_Zero)/2;

    qDebug() << "nL_Zero = " << nL_Zero << "\n";
    qDebug() << "nR_Zero = " << nR_Zero << "\n";
    qDebug() << "N_Zero = " << N_Zero << "\n";


    // Find profile bounds by counting increases/decreases / Nahodim granicy profilya, schitaya uvelicheniya/umensheniya
    // Using second profile (more blurred) for better stability / Ispolzuem vtoroy profil (bolee razmytyy) dlya luchshey stabilnosti
    int decrCounter = 0, incrCounter = 0;
    int XL1 = 0, XR1 = 0;

    // Search right bound / Poisk pravoy granicy
    for(int i = N_Zero + 1; i < n_myu - 1; i++) {
        if(prof2[i] < prof2[i+1]) {
            incrCounter++;
            incrCounter = min(incrCounter, 5);  // Limit to 5 / Ogranichivaem do 5
        }
        if(prof2[i] > prof2[i+1]) {
            decrCounter++;
            decrCounter = min(decrCounter, 5);
        }
        if(decrCounter == 5 && incrCounter == 5) {
            XR1 = min(i, n_myu - 35);;
            break;
        }
    }

    // Search left bound / Poisk levoy granicy
    incrCounter = 0, decrCounter = 0;
    for(int i = N_Zero - 1; i > 1; i--) {
        if(prof2[i] < prof2[i+1]) {
            incrCounter++;
            incrCounter = min(incrCounter, 5);
        }
        if(prof2[i] > prof2[i+1]) {
            decrCounter++;
            decrCounter = min(decrCounter, 5);
        }
        if(decrCounter == 5 && incrCounter == 5) {
            XL1 = max(i, 20);
            break;
        }
    }

    qDebug() << "XL1 = " << XL1 << "\n";
    qDebug() << "XR1 = " << XR1 << "\n";


    // Prepare first profile for interpolation / Podgotavlivaem pervyy profil dlya interpolyacii
    QVector<QPointF> yP1, yP2;
    for(int i = XL1; i < XR1; i++) {
        yP1.push_back({ 1. * i, prof1[i]} );
    }

    int NN = 100;  // New profile length / Novaya dlina profilya

    // Reinterpolate first profile / Pereinterpoliruem pervyy profil
    QVector<QPointF> yyP1 = ImageProcessor::reInterpolateProfile(yP1,NN);
    QVector<double> yy1;
    for(int s1 = 0; s1 < NN; s1++) {
        yy1.push_back(yyP1[s1].y());
    }

    // Show reinterpolated first profile / Pokazyvaem pereinterpolirovannyy pervyy profil
    QVector<double> xProf_yy1, xProf_y1;
    for(int s1 = 0; s1 < NN; s1++) xProf_yy1.push_back(s1);
    for(int s = 0; s < n_myu; s++) {
        xProf_y1.push_back(1.*(s-XL1)/(XR1-XL1)*NN);
    }
    GraphWidget* grwid3 = new GraphWidget();
    grwid3->plotTwoGraphs(xProf_y1, y1, xProf_yy1, yy1);
    grwid3->setWindowTitle(QString("First profile (y1) after reinterpolation (into yy1)"));
    grwid3->show();

    // Calculate scaling factor / Vychislyaem masshtabnyy faktor
    double KF = sqrt((sigma1*sigma1+sigma0*sigma0)/(sigma1*sigma1 + sigma01*sigma01));
    qDebug() << "KF = " << KF << "\n";

    // Search for optimal second profile bounds / Poisk optimalnyh granic vtorogo profilya
    int otstup = 15;  // Search range / Diapozon poiska
    double minraz = INT_MAX;
    int best_XL2, best_XR2;


    for(int XL2 = XL1; XL2 > XL1 - otstup; XL2--) {
        for(int XR2 = XR1; XR2 < XR1 + otstup; XR2++) {

            // Prepare second profile / Podgotavlivaem vtoroy profil
            for(int i = 0; i < n_myu; i++) {
                yP2.push_back({ 1. * i, prof2[i]} );
            }

            // Reinterpolate second profile with current bounds / Pereinterpoliruem vtoroy profil s tekushchimi granicami
            QVector<QPointF> yyP2 = ImageProcessor::reInterpolateProfile(yP2, XL2, XR2, NN);
            QVector<double> yy2;
            for(int s1 = 0; s1 < NN; s1++) {
                yy2.push_back(yyP2[s1].y());
            }

            // Normalize both profiles / Normalizuem oba profilya
            double m1 = 0.0, m2 = 0.0, D1 = 0.0, D2 = 0.0;
            for(int i = 0; i < NN; i++) {
                m1 += yy1[i];
                m2 += yy2[i];
            }
            m1 /= NN + 1;
            m2 /= NN + 1;

            for(int i = 0; i < NN; i++) {
                D1 += (yy1[i] - m1) * (yy1[i] - m1);
                D2 += (yy2[i] - m2) * (yy2[i] - m2);
            }
            D1 = sqrt(D1 / (NN + 1));
            D2 = sqrt(D2 / (NN + 1));

            QVector<double> yyn1, yyn2;
            for(int s1 = 0; s1 < NN; s1++) {
                yyn1.push_back((yy1[s1] - m1) / D1);
                yyn2.push_back((yy2[s1] - m2) / D2);
            }

            // Calculate residual (sum of squared differences) / Vychislyaem nevyazku (summu kvadratov raznostey)
            double MINRAZ = 0.0;
            for(int i = 0; i < NN; i++) {
                MINRAZ += (yyn1[i]-yyn2[i])*(yyn1[i]-yyn2[i]);
            }
            qDebug() << "MINRAZ for (" << XL2 << "," << XR2 << ") = " << MINRAZ << "\n";

            // Keep best bounds / Sohranyaem luchshie granicy
            if(minraz > MINRAZ) {
                minraz = MINRAZ;
                best_XL2 = XL2;
                best_XR2 = XR2;
            }
        }
    }

    qDebug() << "XL2 = " << best_XL2 << " XR2 = " << best_XR2;
    // Final processing with optimal bounds / Finalnaya obrabotka s optimalnymi granicami
    for(int i = 0; i < n_myu; i++) {
        yP2.push_back({ 1. * i, prof2[i]} );
    }

    QVector<QPointF> yyP2 = ImageProcessor::reInterpolateProfile(yP2, best_XL2, best_XR2, NN);
    QVector<double> yy2;
    for(int s1 = 0; s1 < NN; s1++) {
        yy2.push_back(yyP2[s1].y());
    }

    // Final normalization / Finalnaya normalizaciya
    double m1 = 0.0, m2 = 0.0, D1 = 0.0, D2 = 0.0;
    for(int i = 0; i < NN; i++) {
        m1 += yy1[i];
        m2 += yy2[i];
    }
    m1 /= NN + 1;
    m2 /= NN + 1;

    for(int i = 0; i < NN; i++) {
        D1 += (yy1[i] - m1) * (yy1[i] - m1);
        D2 += (yy2[i] - m2) * (yy2[i] - m2);
    }
    D1 = sqrt(D1 / (NN + 1));
    D2 = sqrt(D2 / (NN + 1));

    QVector<double> yyn1, yyn2;
    for(int s1 = 0; s1 < NN; s1++) {
        yyn1.push_back((yy1[s1] - m1) / D1);
        yyn2.push_back((yy2[s1] - m2) / D2);
    }

    // Normalized original profiles / Normalizovannye iskhodnye profili
    QVector<double> yk1, yk2;
    for(int s = 0; s < n_myu; s++) {
        yk1.push_back((y1[s] - m1) / D1);
        yk2.push_back((y2[s] - m2) / D2);
    }

    // Calculate final residual / Vychislyaem finalnuyu nevyazku
    double MINRAZ = 0.0;
    for(int i = 0; i < NN; i++) {
        MINRAZ += (yyn1[i]-yyn2[i])*(yyn1[i]-yyn2[i]);
    }
    qDebug() << "MINRAZ for (" << best_XL2 << "," << best_XR2 << ") = " << MINRAZ << "\n";

    // Calculate shift coefficients / Vychislyaem koefficienty sdviga
    double xa, ya, xb, yb;
    xa = XL1; xb = XR1;
    ya = 1. * best_XL2 - XL1;
    yb = 1. * best_XR2 - XR1;

    qDebug() << "xa = " << xa << "  xb = " << xb;
    qDebug() << "ya = " << ya << "  yb = " << yb << "\n";

    // Find intersection point / Nahodim tochku peresecheniya
    double x00 = (xa*yb - xb*ya) / (yb - ya);
    int n00 = ceil(x00);

    // Convert to mu space / Preobrazuem v prostranstvo mu
    double mu00 = -n_sigma*sigma_myu + x00*(2*n_sigma*sigma_myu)/n_myu;
    qDebug() << "mu00 = " << mu00 << "\n";
    qDebug() << "x00 = " << x00 << "\n";
    qDebug() << "n00 = " << n00 << "\n";

    // Find zero crossings in interpolated profiles / Nahodim perehody cherez nol v interpolirovannyh profilyah
    nR_Zero = N0;
    for(int i = 0; i < NN -  1; i++) {
        if(yy1[i]*yy1[i+1] < 0 && yy1[i] > yy1[i+1]) nL_Zero = i;
        if(yy2[i]*yy2[i+1] < 0 && yy2[i] > yy2[i+1]) nR_Zero = i;
    }
    qDebug() << "nL_Zero = " << nL_Zero << "  nR_Zero = " << nR_Zero << "\n";

    // Final shift and scale coefficients / Finalnye koefficienty sdviga i masshtaba
    double FFF1_1 = mu00;  // Shift coefficient / Koefficient sdviga
    double FFF1_0 = (yb - ya) / (xb - xa);  // Scale coefficient / Koefficient masshtaba
    qDebug() << "FFF1_1 = " << FFF1_1 << "  FFF1_0 = " << FFF1_0 << "\n";

    // Theoretical coefficients for verification / Teoreticheskie koefficienty dlya proverki
    double koef1 = sqrt((sigma1*sigma1 + sigma0*sigma0) / (2*sigma1*sigma1 + sigma0*sigma0));
    double koef2 = 1. * (XR1 - XL1) / (best_XR2 - best_XL2);

    qDebug() << "koef1 = " << koef1 << "  koef2 = " << koef2 << "\n";

    // Transform mu for second profile / Preobrazuem mu dlya vtorogo profilya
    QVector<double> mu_new;
    for(int i = 0; i < n_myu; i++) {
        double val = (mu[i] - FFF1_1) * (1 - FFF1_0) + FFF1_1;
        mu_new.push_back(val);
    }

    // Show aligned profiles / Pokazyvaem vyrovnennye profili
    GraphWidget* grwid7 = new GraphWidget();
    grwid7->plotTwoGraphs(mu, yk1, mu_new, yk2);
    grwid7->setWindowTitle(QString("Profiles yk1 and yk2 after aligning"));
    grwid7->show();

    // Calculate new refined position / Vychislyaem novuyu utochnennuyu poziciyu
    double n_new, m_new;
    n_new = n0 + ex * FFF1_1;
    m_new = m0 + ey * FFF1_1;

    qDebug() << "n = " << n0 << "  m0 = " << m0 << "\n";
    qDebug() << "n_new = " << n_new << "  m_new = " << m_new << "\n";
    qDebug() << "ex = " << ex << "  ey = " << ey << "\n";

    // Create white background image for visualization / Sozdaem izobrazhenie s belym fonom dlya vizualizacii
    QImage tmp = m_dataModel->originalImage();
    QImage whiteIm = tmp;
    for(int i = 0; i < tmp.width(); i++){
        for(int j = 0; j < tmp.height(); j++){
            if(qGray(whiteIm.pixel(i,j)) != 0.) whiteIm.setPixel(i,j,qRgb(255,255,255));
            else whiteIm.setPixel(i,j,qRgb(0,0,0));
        }
    }

    // Mark original and refined points / Otmechaem iskhodnuyu i utochnennuyu tochki
    whiteIm.setPixel(n0, m0, qRgb(0, 0, 255));    // Blue - original point / Siniy - iskhodnaya tochka
    whiteIm.setPixel(n_new, m_new, qRgb(255, 0, 0));    // Red - refined point / Krasnyy - utochnennaya tochka

    // Show final result / Pokazyvaem finalnyy rezultat
    ImageShowcaseWidget* resultWidget = new ImageShowcaseWidget();
    resultWidget->setImage(whiteIm);
    resultWidget->setWindowTitle("Edge Refinement Result");
    resultWidget->show();
}

// ============================================================================
// Image Calculator
// ============================================================================

void MainWindow::on_actionimage_calculator_triggered()
{
    if (curImage.isNull()) {
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    imCalculator->setImage1(curImage, "Current Image");
    imCalculator->show();
    imCalculator->setWindowTitle("Image Calculator");
}

void MainWindow::catch_ImageCalculator(const QImage& image1, const QImage& image2,
                                        QString operation, bool newWindow, bool floatResult)
{
    Q_UNUSED(floatResult);

    if (image1.isNull()) {
        qDebug() << "Error: Image 1 is null";
        return;
    }

    QImage res(image1.width(), image1.height(), QImage::Format_ARGB32);
    Matrix2D<double> ddMat1 = ImageProcessor::fromGrayImage(image1);
    Matrix2D<double> ddMat2 = ImageProcessor::fromGrayImage(image2);

    if (operation == "Add") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::AddWithMaxClamp<double>{}));
    } else if (operation == "Subtract") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::SubtractWithZeroClamp<double>{}));
    } else if (operation == "Multiply") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Multiply<double>{}));
    } else if (operation == "Divide") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Divide<double>{}));
    } else if (operation == "AND") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::BitwiseAND<double>{}));
    } else if (operation == "OR") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::BitwiseOR<double>{}));
    } else if (operation == "XOR") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::BitwiseXOR<double>{}));
    } else if (operation == "Min") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Min<double>{}));
    } else if (operation == "Max") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Max<double>{}));
    } else if (operation == "Average") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Average<double>{}));
    } else if (operation == "Difference") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Difference<double>{}));
    } else if (operation == "Copy") {
        res = image1;
    }

    if (newWindow) {
        ImageShowcaseWidget* showIm = new ImageShowcaseWidget();
        showIm->setImage(res);
        showIm->setWindowTitle("Image Calculator Result");
        showIm->show();
    } else {
        updateImageDisplay();
    }
}

// ============================================================================
// Image Display and Updates
// ============================================================================

void MainWindow::on_imageUpdated(const QImage &newImage)
{
    // Clear matrices
    clearMatrix2D(attMat);
    clearMatrix2D(matForTest);
    clearMatrix2D(imMat);

    curImage = newImage;
    im1 = curImage;
    imMat = ImageProcessor::fromGrayImage(curImage);

    int width = curImage.width();
    int height = curImage.height();

    // Resize matrices
    attMat.resize(width);
    for (int i = 0; i < width; i++) {
        attMat[i].resize(height, 0);
    }

    matForTest.resize(width);
    for (int i = 0; i < width; i++) {
        matForTest[i].resize(height, 0.0);
    }

    updateImageDisplay();
}

void MainWindow::updateImageDisplay()
{
    imageWidget->setImage(curImage);
    imageWidget->show();
}

void MainWindow::showStatisticsDialog(const ImageProcessor::ImageStatistics& stats)
{
    QMessageBox msgBox;
    msgBox.setWindowTitle("Image Statistics");
    msgBox.setText(stats.toString());
    msgBox.setTextFormat(Qt::PlainText);
    msgBox.setStandardButtons(QMessageBox::Ok);

    // Add export buttons
    QPushButton* exportStatsBtn = msgBox.addButton("Export Statistics", QMessageBox::ActionRole);
    QPushButton* exportHistBtn = msgBox.addButton("Export Histogram", QMessageBox::ActionRole);

    msgBox.exec();

    if (msgBox.clickedButton() == exportStatsBtn) {
        on_actionExportStatistics_triggered();
    } else if (msgBox.clickedButton() == exportHistBtn) {
        on_actionExportHistogram_triggered();
    }
}

void MainWindow::showProfileDialog(const QVector<double>& profile, const QString& title)
{
    if (profile.isEmpty()) return;

    QDialog dialog(this);
    dialog.setWindowTitle(title);
    dialog.resize(600, 400);

    QVBoxLayout* layout = new QVBoxLayout(&dialog);

    GraphWidget* graph = new GraphWidget(&dialog);
    QVector<double> x(profile.size());
    for (int i = 0; i < profile.size(); i++) x[i] = i;
    graph->plotGraph(x, profile);
    graph->setAxisLabels("Distance (pixels)", "Intensity");
    graph->setTitle(title);

    layout->addWidget(graph);

    QPushButton* closeBtn = new QPushButton("Close", &dialog);
    layout->addWidget(closeBtn);

    connect(closeBtn, &QPushButton::clicked, &dialog, &QDialog::accept);

    dialog.exec();
}

// ============================================================================
// Unused Actions (Placeholders)
// ============================================================================

void MainWindow::on_actionpcr_edge_detection_triggered()
{
    QMessageBox::information(this, "PCR Edge Detection",
        "PCR edge detection is not yet implemented. Use test_001 for refinement.");
}
