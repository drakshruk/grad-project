#include "mainwindow.h"
#include "ui_mainwindow.h"

QPoint mPos;
QVector<QPoint> selectedEdge;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_dataModel(new AppDataModel(this))
{
    ui->setupUi(this);
    imCalculator = new ImageCalculator();
    imageWidget = new ImageShowcaseWidget();

    connect(m_dataModel, &AppDataModel::currentImageChanged, this, &MainWindow::on_imageUpdated);
    connect(imCalculator,&ImageCalculator::throw_imageCalculator , this, &MainWindow::catch_ImageCalculator);
    connect(imageWidget, &ImageShowcaseWidget::imageClicked, this, &MainWindow::on_showcaseWidget_clicked);

    m_dataModel->setBlueMax(255);
    m_dataModel->setRedMax(255);
    m_dataModel->setRadius(10);
    m_dataModel->setSigma(5);

    getDataDialog = new Dialog();
    this->setFixedSize(550,50);
    GrWid = new GraphWidget();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actiongaussian_blur_triggered()
{
    this->getDataDialog->setWindowTitle("Write sigma coefficient value");
    this->getDataDialog->setPlaceholderText("Sigma value");
    if(getDataDialog->exec())
    {
        m_dataModel->setSigma(getDataDialog->value);
        m_dataModel->setRadius(m_dataModel->sigma());
        curImage = ImageProcessor::convImage(curImage, ImageProcessor::getGauss(m_dataModel->radius(),m_dataModel->radius(),m_dataModel->sigma()));

        imageWidget->setImage(curImage);
        imageWidget->show();
        im1 = curImage;
    }
}

void MainWindow::on_actiongaussian_edge_detection_triggered()
{
    this->getDataDialog->setWindowTitle("Write sigma coefficient value");
    this->getDataDialog->setPlaceholderText("Sigma value");
    getDataDialog->exec();
    m_dataModel->setSigma(getDataDialog->value);
    m_dataModel->setRadius(3*m_dataModel->sigma());
    curImage = ImageProcessor::gaussianEdgeDetection(ImageProcessor::fromGrayImage(curImage), m_dataModel->sigma(), m_dataModel->radius(), attMat, matForTest);

    edgeSelectionMode = true;

    imageWidget->setImage(curImage);
    imageWidget->show();

    edgeSelectionMode = true;
}

void MainWindow::on_actionLaplacian_edge_detection_triggered()
{
    this->getDataDialog->setWindowTitle("Write sigma coefficient value");
    this->getDataDialog->setPlaceholderText("Sigma value");
    getDataDialog->exec();
    m_dataModel->setSigma(getDataDialog->value);
    m_dataModel->setRadius(3*m_dataModel->sigma());

    Matrix2D<double> mat1;
    Matrix2D<double> mat2;

    mat1 = ImageProcessor::fromGrayImage(curImage);
    mat2 = ImageProcessor::convMat(mat1, ImageProcessor::getLapl(m_dataModel->radius(), m_dataModel->radius(), m_dataModel->sigma()));
    matForTest = mat2;
    mat1 = ImageProcessor::elementWiseOperation(mat1,mat2,MatrixLambdas::Subtract<double>{});
    matForTest = mat1;
    mat1 = ImageProcessor::findEdges(mat1, 3);
    curImage = ImageProcessor::toGrayImage(mat1);
    edgeSelectionMode = true;
}

void MainWindow::on_imageUpdated(const QImage &newImage)
{
    clearMatrix2D(attMat);
    clearMatrix2D(matForTest);
    clearMatrix2D(imMat);

    curImage = newImage;
    im1 = curImage;

    imMat = ImageProcessor::fromGrayImage(curImage);

    int width = curImage.width();
    int height = curImage.height();

    attMat.resize(width);
    for(int i = 0; i < width; i++) {
        attMat[i].resize(height, 0);
    }

    matForTest.resize(width);
    for(int i = 0; i < width; i++) {
        matForTest[i].resize(height, 0.0);
    }

    imageWidget->setImage(curImage);
    imageWidget->show();
}

void MainWindow::on_actionpcr_edge_detection_triggered()
{

}

void MainWindow::on_actionopen_file_triggered()
{
    QFileDialog* openImDialog = new QFileDialog(this);
    openImDialog->setFileMode(QFileDialog::AnyFile);
    openImDialog->setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));

    QImage newIm = QImage(openImDialog->getOpenFileName());
    m_dataModel->setCurrentImage(newIm);
    m_dataModel->setOriginalImage(newIm);
    this->show();
}

void MainWindow::on_actionsave_file_triggered()
{
    QFileDialog* saveImDialog = new QFileDialog(this);
    saveImDialog->setFileMode(QFileDialog::AnyFile);
    saveImDialog->setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));
    QString saveFileName = saveImDialog->getSaveFileName();
    curImage.save(saveFileName);
}

void MainWindow::on_showcaseWidget_clicked(const QPoint &imagePosition)
{
    mPos.setX(imagePosition.x());
    mPos.setY(imagePosition.y());
    if(edgeSelectionMode){
        int cnt = 0;
        for(int i = 0; i < curImage.width(); i++){
            for(int j = 0; j < curImage.height(); j++){
                if(attMat[i][j] == int(attribute::isSelectedEdge)){
                    attMat[i][j] = int(attribute::isEdge);
                    curImage.setPixel(i,j,qRgb(255,255,255));
                }
                if(attMat[i][j] == int(attribute::isEdge)){
                    cnt++;
                }
            }
        }
        selectedEdge.clear();
        ImageProcessor::selectEdge(mPos, attMat, selectedEdge);
        qDebug() << "selected edge size = " << selectedEdge.size();
        qDebug() << "counter = " << cnt;
        m_dataModel->setSelectedEdge(selectedEdge);
        for(int i = 0; i < selectedEdge.size(); i++){
            curImage.setPixel(selectedEdge[i],qRgb(0,255,255));
        }
        imageWidget->setImage(curImage);
        imageWidget->show();
    }
}

void MainWindow::catch_ImageCalculator(const QImage& image1, const QImage& image2, QString operation, bool newWindow, bool floatResult)
{
    Q_UNUSED(floatResult);
    QImage res(image1.width(),image1.height(),QImage::Format_ARGB32);
    Matrix2D<double> ddMat1 = ImageProcessor::fromGrayImage(image1);
    Matrix2D<double> ddMat2 = ImageProcessor::fromGrayImage(image2);
    if(operation == "Add"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::AddWithMaxClamp<double>{}));
    } else if(operation == "Substract"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::SubtractWithZeroClamp<double>{}));
    } else if(operation == "Multiply"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::Multiply<double>{}));
    } else if(operation == "Divide"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::Divide<double>{}));
    } else if(operation == "AND"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::BitwiseAND<double>{}));
    } else if(operation == "OR"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::BitwiseOR<double>{}));
    } else if(operation == "XOR"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::BitwiseXOR<double>{}));
    } else if(operation == "Min"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::Min<double>{}));
    } else if(operation == "Max"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::Max<double>{}));
    } else if(operation == "Average"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::Average<double>{}));
    } else if(operation == "Difference"){
        res = ImageProcessor::toGrayImage(ImageProcessor::elementWiseOperation(ddMat1,ddMat2,MatrixLambdas::Difference<double>{}));
    } else if(operation == "Copy"){
        res = image1;
    } else if(operation == "Transparent-zero"){

    }
    if(newWindow) {
        ImageShowcaseWidget* showIm = new ImageShowcaseWidget();
        showIm->setImage(res);
        showIm->show();
    }
    else imageWidget->setImage(res);
    clearMatrix2D(ddMat1);
    clearMatrix2D(ddMat2);
}

void MainWindow::on_actiongradient_X_and_Y_triggered()
{
    Matrix2D<double> xGrad, yGrad;
    xGrad = ImageProcessor::fromGrayImage(curImage);
    xGrad = ImageProcessor::convMat(xGrad, ImageProcessor::getXGradCore(iRad,iRad,dSigma));
    ImageProcessor::toBlueRedImage(xGrad, 255., 255.).save("D:/projects/ERWS/Grad/XGradImage.png");
    yGrad = ImageProcessor::fromGrayImage(curImage);
    yGrad = ImageProcessor::convMat(yGrad, ImageProcessor::getYGradCore(iRad,iRad,dSigma));
    ImageProcessor::toBlueRedImage(yGrad, 255., 255.).save("D:/projects/ERWS/Grad/YGradImage.png");
}

void MainWindow::on_actiondraw_profile_triggered()
{
    QVector<double> x, y;
    for(int i = 0; i < curImage.height(); i++){
        x.push_back(i);
        y.push_back(qGray(curImage.pixel(curImage.width()/2,i)));
    }
    GrWid->plotGraph(x,y);
    GrWid->show();
}

void MainWindow::on_actionTwo_hollows_triggered()
{
    QImage setImage = ImageProcessor::sampleTwoHollows();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);
}

void MainWindow::on_actionTwo_hollows_big_triggered()
{
    QImage setImage = ImageProcessor::sampleTwoHollowsBig();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);
}

void MainWindow::on_actionimage_calculator_triggered()
{
    imCalculator->show();
    imCalculator->setWindowTitle("Image calculator");
}

void MainWindow::on_actionsharpen_triggered()
{
    Matrix2D<double> kernel = {{-1,-1,-1},
                               {-1, 9,-1},
                               {-1,-1,-1}};
    curImage = ImageProcessor::convImage(curImage, kernel);

    imageWidget->setImage(curImage);
    imageWidget->show();
    im1 = curImage;
}

void MainWindow::on_actiontest_triggered()
{
    // Check if edge selection mode is active / Proveryaem aktivirovan li rezhim vybora granic
    if(!edgeSelectionMode) {
        qDebug() << "Error: Edge selection mode is off";
        return;
    }

    // Collect all selected edge points / Sobiraem vse vybrannye granichnye tochki
    QVector<QPoint> selectedEdgePoints;
    for(int i = 0; i < (int)attMat.size(); i++) {
        for(int j = 0; j < (int)attMat[i].size(); j++) {
            if(attMat[i][j] == static_cast<int>(attribute::isSelectedEdge)) {
                selectedEdgePoints.append(QPoint(i, j));
            }
        }
    }

    // Check if any points were selected / Proveryaem, byli li vybrany tochki
    if(selectedEdgePoints.isEmpty()) {
        qDebug() << "No edge points selected";
        return;
    }

    qDebug() << "Refining" << selectedEdgePoints.size() << "edge points";

    // Create white background image for visualization / Sozdaem izobrazhenie s belym fonom dlya vizualizacii
    QImage whiteImage = m_dataModel->originalImage();
    for(int i = 0; i < whiteImage.width(); i++) {
        for(int j = 0; j < whiteImage.height(); j++) {
            if(qGray(whiteImage.pixel(i,j)) != 0) {
                whiteImage.setPixel(i, j, qRgb(255, 255, 255));
            } else {
                whiteImage.setPixel(i, j, qRgb(0, 0, 0));
            }
        }
    }

    // Image dimensions / Razmery izobrazheniya
    int NX = m_dataModel->originalImage().width();
    int NY = m_dataModel->originalImage().height();

    // Gaussian blur parameters / Parametry razmytiya Gaussa
    double sigma0 = 8.0, sigma01 = sqrt(sigma0*sigma0 + 4.*4.);
    int NS = 2*int(4*sigma01+0.5);

    // Two different sigma values for edge detection / Dva raznyh znacheniya sigma dlya detektirovaniya granic
    double sigma1 = 4., sigma2 = 8.;
    double sigma_max = max(sigma1, sigma2);

    // Prepare matrices (once for all points) / Podgotavlivaem matricy (odin raz dlya vseh tochek)
    Matrix2D<double> A0 = ImageProcessor::fromGrayImage(m_dataModel->originalImage());
    Matrix2D<double> B1 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma0));
    Matrix2D<double> B01 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma01));
    Matrix2D<double> A = B1;

    // Update kernel size for max sigma / Obnovlyaem razmer yadra dlya maksimalnoy sigma
    NS = 2*int(4*sigma_max+0.5);

    // Calculate gradients / Vychislyaem gradienty
    Matrix2D<double> B = A;
    Matrix2D<double> GxW1 = ImageProcessor::convMat(B, ImageProcessor::getXGradCore(NS, NS, sigma1));
    Matrix2D<double> GyW1 = ImageProcessor::convMat(B, ImageProcessor::getYGradCore(NS, NS, sigma1));

    // Profile parameters / Parametry profilya
    int n_sigma = 10;           // Number of sigma steps / Kolichestvo shagov sigma
    int n_myu = 200;             // Number of mu points / Kolichestvo tochek mu
    double sigma_myu = min(sigma1, sigma2);  // Sigma for profile / Sigma dlya profilya
    int NN = 100;
    int otstup = 15;

    // Generate mu values (same for all points) / Generiruem znacheniya mu (odinakovye dlya vseh tochek)
    QVector<double> mu;
    for(int s = 0; s < n_myu; s++) {
        double val = -n_sigma*sigma_myu + s*(2*n_sigma*sigma_myu)/n_myu;
        mu.push_back(val);
    }

    // Vector to store refined positions / Vektor dlya hraneniya utochnennyh poziciy
    QVector<QPointF> refinedPositions;

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

    // Process each selected point / Obrabatyvaem kazhduyu vybrannuyu tochku
    for(int pointIdx = 0; pointIdx < selectedEdgePoints.size(); pointIdx++) {
        int n0 = selectedEdgePoints[pointIdx].x(), m0 = selectedEdgePoints[pointIdx].y();

        double ex = GxW1[n0][m0];
        double ey = GyW1[n0][m0];
        double gradMag = sqrt(ex*ex+ ey*ey);
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

        // Visualization on white image / Vizualizaciya na belom izobrazhenii
        // Original point - blue / Iskhodnaya tochka - sinii
        if(n0 >= 0 && n0 < whiteImage.width() && m0 >= 0 && m0 < whiteImage.height()) {
            whiteImage.setPixel(n0, m0, qRgb(0, 0, 255));
        }

        // Refined point - red / Utochnennaya tochka - krasnyi
        if(n_new >= 0 && n_new < whiteImage.width() && m_new >= 0 && m_new < whiteImage.height()) {
            whiteImage.setPixel(n_new, m_new, qRgb(255, 0, 0));
        }

        qDebug() << "Point" << pointIdx + 1 << " / " << selectedEdgePoints.size()
                 << "refined: (" << n0 << "," << m0 << ") -> ("
                 << n_new << "," << m_new << ") shift =" << result.FFF1_1;
    }

    // Statistics for all points / Statistika po vsem tochkam
    qDebug() << "\n=== REFINEMENT SUMMARY / STATISTIKA UTOCHNENIYA ===";
    qDebug() << "Total points processed / Vsego obrabotano tochek:" << refinedPositions.size();

    double avgShiftX = 0.0, avgShiftY = 0.0;
    for(int i = 0; i < selectedEdgePoints.size(); i++) {
        double shiftX = refinedPositions[i].x() - selectedEdgePoints[i].x();
        double shiftY = refinedPositions[i].y() - selectedEdgePoints[i].y();
        avgShiftX += shiftX;
        avgShiftY += shiftY;
        qDebug() << "Point" << i+1 << ": original / iskhodnaya (" << selectedEdgePoints[i].x()
                 << "," << selectedEdgePoints[i].y() << ") -> refined / utochnennaya ("
                 << refinedPositions[i].x() << "," << refinedPositions[i].y()
                 << ") shift / sdvig (" << shiftX << "," << shiftY << ")";
    }
    avgShiftX /= refinedPositions.size();
    avgShiftY /= refinedPositions.size();
    qDebug() << "Average shift / Srednii sdvig: (" << avgShiftX << "," << avgShiftY << ")";

    // Show final image / Pokazyvaem finalnoe izobrazhenie
    ImageShowcaseWidget* resultWidget = new ImageShowcaseWidget();
    resultWidget->setImage(whiteImage);
    resultWidget->setWindowTitle(QString("Edge Refinement Result / Rezultat utochneniya granic - %1 points").arg(refinedPositions.size()));
    resultWidget->show();

    // Save result / Sohranyaem rezultat
    resultWidget->getImage().save("testing.png");

    qDebug() << "Edge refinement completed / Utochnenie granic zaversheno";
}

void MainWindow::on_actiontest_002_triggered()
{
    // Generate test image with two hollows / Sozdaem testovoe izobrazhenie s dvumya polostyami
    QImage setImage = ImageProcessor::sampleTwoHollows();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);

    // Image dimensions / Razmery izobrazheniya
    int NX = setImage.width();
    int NY = setImage.height();

    // Gaussian blur parameters / Parametry razmytiya Gaussa
    double sigma0 = 8.0, sigma01 = sqrt(sigma0*sigma0 + 4.*4.);
    int NS = 2*int(4*sigma01+0.5);
    qDebug() << "NS = " << NS << "\n";

    // Two different sigma values for edge detection / Dva raznyh znacheniya sigma dlya detektirovaniya granic
    double sigma1 = 4., sigma2 = 8.;
    double sigma_max = max(sigma1, sigma2);

    // Convert images to matrices and apply Gaussian blur / Preobrazuem izobrazheniya v matricy i primenyaem razmytie Gaussa
    Matrix2D<double> A0 = ImageProcessor::fromGrayImage(m_dataModel->currentImage());
    Matrix2D<double> B1 = ImageProcessor::convMat(A0,ImageProcessor::getGauss(NS,NS,sigma0));
    Matrix2D<double> B01 = ImageProcessor::convMat(A0,ImageProcessor::getGauss(NS,NS,sigma01));
    Matrix2D<double> A = B1;

    // Statistics of blurred image / Statistika razmytogo izobrazheniya
    double A_max = max(A);
    double A_min= min(A);
    int NSIZE = A.size();

    qDebug() << "max(B1) = " << A_max << "\n";
    qDebug() << "min(B1) = " << A_min << "\n";

    // Show profiles before and after convolution / Pokazyvaem profili do i posle svertki
    GraphWidget* grwid1 = new GraphWidget();
    QVector<double> profileA, profileA0;
    QVector<double> profileX;
    for(int i = 0; i < NSIZE; i++) {
        profileA.push_back(A[100][i]);
        profileA0.push_back(A0[100][i]);
        profileX.push_back(i);
    }
    grwid1->plotTwoGraphs(profileX, profileA, profileX, profileA0);
    grwid1->setWindowTitle(QString("Profiles before and after convolution"));
    grwid1->show();

    // Update kernel size for max sigma / Obnovlyaem razmer yadra dlya maksimalnoy sigma
    NS = 2*int(4*sigma_max+0.5);
    qDebug() << "NS = " << NS << "\n";

    // Laplacian kernels for edge detection / Yadra Laplasa dlya detektirovaniya granic
    Matrix2D<double> MS1 = ImageProcessor::getLapl(NS,NS,sigma1);
    Matrix2D<double> MS2 = ImageProcessor::getLapl(NS,NS,sigma2);

    // Perform Gaussian edge detection / Vypolnyaem detektirovanie granic po Gaussu
    curImage = ImageProcessor::gaussianEdgeDetection(ImageProcessor::fromGrayImage(curImage), sigma_max, NS, attMat, matForTest);
    edgeSelectionMode = true;
    imageWidget->setImage(curImage);
    imageWidget->show();

    // Calculate gradients / Vychislyaem gradienty
    Matrix2D<double> B = A;
    Matrix2D<double> GxW1 = ImageProcessor::convMat(B, ImageProcessor::getXGradCore(NS, NS, sigma1));
    Matrix2D<double> GyW1 = ImageProcessor::convMat(B, ImageProcessor::getYGradCore(NS, NS, sigma1));

    // Gradient statistics / Statistika gradientov
    double X_max = max(GxW1), X_min = min(GxW1);
    double Y_max = max(GyW1), Y_min = min(GyW1);
    qDebug() << "x grad min : " << X_min << ", max : " << X_max << ";\n";
    qDebug() << "y grad min : " << Y_min << ", max : " << Y_max << ";\n";

    //****************************************************************************************************
    // Select test point and normalize its gradient / Vyberaem testovuyu tochku i normalizuem ee gradient
    // Start refining position / Nachynaem utochnenie positsii
    //****************************************************************************************************
//    int n0 = 109, m0 = 87;
//    int n0 = 102, m0 = 89;
//    int n0 = 119, m0 = 86;
//    int n0 = 116, m0 = 84;
//    int n0 = 111, m0 = 86;
//    int n0 = 90, m0 = 87;
    int n0 = 78, m0 = 83;
    double ex = GxW1[n0][m0];
    double ey = GyW1[n0][m0];
    double gradMag = sqrt(ex*ex+ ey*ey);
    if (gradMag > 1e-6) {
        ex /= gradMag;
        ey /= gradMag;
    } else {
        ex = 1.0;
        ey = 0.0;
    }
    qDebug() << "x grad at (n0,m0): " << ex << ";\n";
    qDebug() << "y grad at (n0,m0): " << ey << ";\n";

    // Calculate gradient magnitude map / Vychislyaem kartu velichiny gradienta
    Matrix2D<double> MGW1(NSIZE,std::vector<double>(NSIZE,0.0));
    for(int i = 0; i < NSIZE; i++) {
        for(int j = 0; j < NSIZE; j++) {
            MGW1[i][j] = sqrt(GxW1[i][j]*GxW1[i][j] + GyW1[i][j]*GyW1[i][j]);
        }
    }

    // Show gradient magnitude image / Pokazyvaem izobrazhenie velichiny gradienta
    QImage mgw_image = ImageProcessor::toGrayImage(MGW1);
    ImageShowcaseWidget* mgw_widget = new ImageShowcaseWidget();
    mgw_widget->setImage(mgw_image);
    mgw_widget->show();
    qDebug() << "mgw image at (n0,m0): " << MGW1[n0][m0] << ";\n";
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
