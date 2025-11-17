#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "utility.cpp"

/*
 *
 */
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

    m_dataModel->setBlueMax(255);
    m_dataModel->setRedMax(255);
    m_dataModel->setRadius(10);
    m_dataModel->setSigma(5);

    getDataDialog = new Dialog();
//    curImage = QImage("D:/projects/ERWS/Grad/files/fto2.png");
//    m_dataModel->setCurrentImage(curImage);
    this->setFixedSize(550,50);
//    ui->label->setPixmap(pixmap->fromImage(curImage));
//    im1 = curImage;
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
//        ui->label->setPixmap(pixmap->fromImage(curImage));

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

//    m_dataModel->setCurrentImage(curImage);
    ui->label->setPixmap(pixmap->fromImage(curImage));

    imageWidget->setImage(curImage);
    imageWidget->show();
//    curImage.save("D:/projects/ERWS/Grad/edgeImage.png");
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
    ui->label->setPixmap(pixmap->fromImage(curImage));
    curImage.save("D:/projects/ERWS/Grad/edgeImage.png");
    edgeSelectionMode = true;
}

void MainWindow::on_imageUpdated(const QImage &newImage)
{
    clearMatrix2D(attMat);
    clearMatrix2D(matForTest);
    clearMatrix2D(imMat);

    curImage = newImage;
    im1 = curImage;
    this->setFixedSize(newImage.width(),newImage.height()+100);
    ui->label->setPixmap(pixmap->fromImage(newImage));

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
    /*
     * median and dispersion
     * fix kernels, sum and divide
     *
     */
    if(edgeSelectionMode) { edgeSelectionMode = false; qDebug() << "edge selection off, starting edge refining"; }

    int nSig = 8;
    this->getDataDialog->setWindowTitle("Write the amount of sigms in the profile");
    this->getDataDialog->setPlaceholderText("Amount of sigms");
    if(getDataDialog->exec()){
        nSig = int(getDataDialog->value);
    }
    qDebug() << "1.1";

    Matrix2D<double> originalMat = ImageProcessor::fromGrayImage(m_dataModel->originalImage());

    ImageShowcaseWidget* im = new ImageShowcaseWidget();
    im->setImage(m_dataModel->originalImage());
    im->show();

    qDebug() << "1.11";
    double sigma1 = m_dataModel->sigma();
    double sigma2 = sigma1*1.5;
    qDebug() << "1.12";
    Matrix2D<double> blurredMat1 = ImageProcessor::convMat(originalMat,ImageProcessor::getGauss(m_dataModel->radius(), m_dataModel->radius(), sigma1));
    qDebug() << "1.13";
    Matrix2D<double> blurredMat2 = ImageProcessor::convMat(originalMat,ImageProcessor::getGauss(m_dataModel->radius(), m_dataModel->radius(), sigma2));
    qDebug() << "1.2";

    Matrix2D<double> xGradMat = ImageProcessor::convMat(originalMat, ImageProcessor::getXGradCore(m_dataModel->radius(), m_dataModel->radius(), sigma1));
    Matrix2D<double> yGradMat = ImageProcessor::convMat(originalMat, ImageProcessor::getYGradCore(m_dataModel->radius(), m_dataModel->radius(), sigma1));

    ImageShowcaseWidget* im1 = new ImageShowcaseWidget();
    im1->setImage(ImageProcessor::toBlueRedImage(xGradMat,255,255));
    im1->show();

    ImageShowcaseWidget* im2 = new ImageShowcaseWidget();
    im2->setImage(ImageProcessor::toBlueRedImage(yGradMat,255,255));
    im2->show();

    ImageShowcaseWidget* im3 = new ImageShowcaseWidget();
    im3->setImage(ImageProcessor::toGrayImage(blurredMat1));
    im3->show();

    ImageShowcaseWidget* im4 = new ImageShowcaseWidget();
    im4->setImage(ImageProcessor::toGrayImage(blurredMat2));
    im4->show();

    blurredMat1 = ImageProcessor::fromGrayImage(im3->getImage().toImage());
    blurredMat2 = ImageProcessor::fromGrayImage(im4->getImage().toImage());
    qDebug() << "1.3";

    QImage tmp = im->getImage().toImage();

    QImage whiteIm = tmp;
    for(int i = 0; i < tmp.width(); i++){
        for(int j = 0; j < tmp.height(); j++){
            if(qGray(whiteIm.pixel(i,j)) != 0.) whiteIm.setPixel(i,j,qRgb(255,255,255));
            else whiteIm.setPixel(i,j,qRgb(0,0,0));
        }
    }

    GraphWidget* wid1 = new GraphWidget();
    GraphWidget* wid2 = new GraphWidget();
    GraphWidget* debugWid = new GraphWidget();

    int searchRange = 15;
    if(getDataDialog->exec()){
        searchRange = int(getDataDialog->value);
    }

    for(int i = 9; i < selectedEdge.size(); i++)
    {

        double xGrad = 0., yGrad = 0.;
        xGrad = xGradMat[selectedEdge[i].x()][selectedEdge[i].y()];
        yGrad = yGradMat[selectedEdge[i].x()][selectedEdge[i].y()];
//        qDebug() << " xGrad: " << xGrad << " yGrad: " << yGrad;
        double gradMag = sqrt(xGrad*xGrad + yGrad*yGrad);
        if (gradMag > 1e-6) {
            xGrad /= gradMag;
            yGrad /= gradMag;
        } else {
            xGrad = 1.0;
            yGrad = 0.0;
        }

//        qDebug() << "2.1";

        ProfileParameters params;


//        qDebug() << "2.2";
//        qDebug() << " xGrad: " << xGrad << " yGrad: " << yGrad;
        params.numSigma = nSig;
        params.sigma = sigma1;
        params.x0 = selectedEdge[i].x();
        params.y0 = selectedEdge[i].y();
        params.xGradient = xGrad;
        params.yGradient = yGrad;

        ProfileResult profile1 = ImageProcessor::buildProfile002(params, blurredMat1, 1);
        ProfileResult profile2 = ImageProcessor::buildProfile002(params, blurredMat2, 1);
        if(i <= 2){

            QVector<double> x, y;
            for(int j = 0; j < profile1.points.size(); j++){
                x.push_back(profile1.points[j].x());
                y.push_back(profile1.points[j].y());
            }
            wid1->plotGraph(x,y); wid1->show();

            x.clear(); y.clear();
            for(int j = 0; j < profile2.points.size(); j++){
                x.push_back(profile2.points[j].x());
                y.push_back(profile2.points[j].y());
            }
            wid2->plotGraph(x,y); wid2->show();

        }
//        qDebug() << "2.3";

        std::pair<int,int> bounds = ImageProcessor::findProfileBounds(profile1.points, 5);
        int left1 = bounds.first, right1 = bounds.second;

//        qDebug() << "2.4";

        QVector<QPointF> referenceProfile = profile1.points;
        referenceProfile = ImageProcessor::cutProfile(referenceProfile, left1, right1);
        referenceProfile = ImageProcessor::reInterpolateProfile(profile1.points, nSig*dSigma);

//        qDebug() << "2.5";


        std::pair<int,int> bounds2 = ImageProcessor::alignProfiles(referenceProfile, profile2.points, left1, right1, searchRange);
        int left2 = bounds2.first, right2 = bounds2.second;

        if(i == 1){
            QVector<QPointF> debugProf = ImageProcessor::cutProfile(profile2.points, left2, right2);
            debugProf = ImageProcessor::reInterpolateProfile(debugProf, referenceProfile.size());
            QVector<double> x, y;
            for(int j = 0; j < debugProf.size(); j++){
                x.push_back(debugProf[j].x());
                y.push_back(debugProf[j].y());
            }
            debugWid->plotGraph(x,y);
            debugWid->show();
        }
//        qDebug() << "2.6";

        double leftShift = abs(left1 - left2);
        double rightShift = abs(right1 - right2);

//        qDebug() << " profile1 size: " << profile1.points.size() << " profile2 size: " << profile2.points.size();
//        qDebug() << " left1: " << left1 << " right1: " << right1;
//        qDebug() << " left2: " << left2 << " right2: " << right2;
//        qDebug() << " left shift: " << leftShift << " right shift: " << rightShift;
//        double x = (left1*rightShift - right1*leftShift)/(leftShift - rightShift);

//        double mu = nSig*dSigma*(-1. + 2.*x/nSig/dSigma);

        double shift = (leftShift + rightShift) / 2.0;

//        int xNewOriginal = static_cast<int>(mu+0.5);
//        int yNewOriginal = (rightShift-leftShift)/(right1-left1);

        int xNew = selectedEdge[i].x() + shift*xGrad;
        int yNew = selectedEdge[i].y() + shift*yGrad;

        qDebug() << " New x: " << xNew << " New y: " << yNew;
        qDebug() << " Old x: " << selectedEdge[i].x() << " Old y: " << selectedEdge[i].y();
//        qDebug() << " new x original: " << xNewOriginal << " new y original: " << yNewOriginal;
        qDebug() << " residual: " << ImageProcessor::calculateResidualOfProfiles(
                        ImageProcessor::reInterpolateProfile(ImageProcessor::cutProfile(profile1.points,left1,right1),100),
                        ImageProcessor::reInterpolateProfile(ImageProcessor::cutProfile(profile2.points,left2,right2),100));
        tmp.setPixel(xNew, yNew, qRgb(255,0,0));
        tmp.setPixel(selectedEdge[i].x(), selectedEdge[i].y(), qRgb(0,0,255));
        whiteIm.setPixel(xNew, yNew, qRgb(255,0,0));
        whiteIm.setPixel(selectedEdge[i].x(), selectedEdge[i].y(), qRgb(0,0,255));
        im->setImage(tmp);
    }

    ImageShowcaseWidget* im5 = new ImageShowcaseWidget();
    im5->setImage(whiteIm);
    im5->show();

    m_dataModel->setCurrentImage(im->getImage().toImage());
    qDebug() << "edge size: " << selectedEdge.size();
}

void MainWindow::on_actionopen_file_triggered()
{
    QFileDialog* openImDialog = new QFileDialog(this);
    openImDialog->setFileMode(QFileDialog::AnyFile);
    openImDialog->setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));

    QImage newIm = QImage(openImDialog->getOpenFileName());
    m_dataModel->setCurrentImage(newIm);
//    setImage(QImage(openImDialog->getOpenFileName()));
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

void MainWindow::mousePressEvent(QMouseEvent *e)
{
    mPos.setX(e->pos().x()-10);
    mPos.setY(e->pos().y()-55);
    qDebug() << mPos.x() << " " << mPos.y() << " " << attMat[mPos.x()][mPos.y()];
    qDebug() << qGray(curImage.pixel(mPos));
    if(edgeSelectionMode){
        for(int i = 0; i < curImage.width(); i++){
            for(int j = 0; j < curImage.height(); j++){
                if(attMat[i][j] == int(attribute::isSelectedEdge)){
                    attMat[i][j] = int(attribute::isEdge);
                    curImage.setPixel(i,j,qRgb(255,255,255));
                }
            }
        }
        selectedEdge.clear();
        ImageProcessor::selectEdge(mPos, attMat, selectedEdge);
        m_dataModel->setSelectedEdge(selectedEdge);
        for(int i = 0; i < selectedEdge.size(); i++){
            curImage.setPixel(selectedEdge[i],qRgb(0,255,255));
        }
        imageWidget->setImage(curImage);
        imageWidget->show();
//        ui->label->setPixmap(pixmap->fromImage(curImage));
    }
}

void MainWindow::catch_ImageCalculator(const QImage& image1, const QImage& image2, QString operation, bool newWindow, bool floatResult)
{
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

void MainWindow::on_actiontest_triggered()
{
    Matrix2D<double> mat1;
    mat1 = ImageProcessor::fromGrayImage(curImage);
    Matrix2D<double> matSig1, matSig2;
    Matrix2D<double> xGrad, yGrad;
    qDebug() << 1.1;
    xGrad = matForTest;
    qDebug() << 1.2;
    xGrad = ImageProcessor::convMat(xGrad, ImageProcessor::getXGradCore(iRad,iRad,dSigma));
    qDebug() << 1.3;
    yGrad = matForTest;
    qDebug() << 1.4;
    yGrad = ImageProcessor::convMat(yGrad, ImageProcessor::getYGradCore(iRad,iRad,dSigma));
    qDebug() << 1.5;
    double xGradient = 0., yGradient = 0.;

    for(int i = 0; i < curImage.width(); i++)
    {
        for(int j = 0; j < curImage.height(); j++)
        {
            if(attMat[i][j] == static_cast<int>(attribute::isEdge)) mat1[i][j] = 255.;
            else mat1[i][j] = 0.;
        }
    }

    QImage imSet = ImageProcessor::toBlueRedImage(matForTest, 255, 255);
    imageWidget->setImage(curImage);
//    ui->label->setPixmap(pixmap->fromImage(imSet));
    curImage = imSet;
    xGradient = xGrad[mPos.x()][mPos.y()];
    yGradient = yGrad[mPos.x()][mPos.y()];
    double tmp = sqrt(xGradient*xGradient + yGradient*yGradient);
    xGradient /= tmp; yGradient /= tmp;
    int nSig = 8;
    this->getDataDialog->setWindowTitle("Write the amount of sigms in the profile");
    this->getDataDialog->setPlaceholderText("Amount of sigms");
    if(getDataDialog->exec()){
        qDebug() << 1;
        nSig = int(getDataDialog->value);
    }
    qDebug() << 2;
    ProfileParameters params{mPos.x(), mPos.y(), dSigma, xGradient, yGradient, nSig};
/*
    params.numSigma = nSig;
    params.sigma = dSigma;
    params.x0 = mPos.x();
    params.y0 = mPos.y();
    params.xGradient = xGradient;
    params.yGradient = yGradient;
*/
    ProfileResult profile = ImageProcessor::buildProfile002(params, matForTest, 3);
    qDebug() << "xGrad" << xGradient << "yGrad" << yGradient;
    qDebug() << 3 << "nSig" << nSig;
    qDebug() << profile.indices.size() << profile.points.size() << "===============";
    for(int i = 0; i < profile.indices.size(); i++){
         qDebug() << profile.indices[i].x() << profile.indices[i].y();
        curImage.setPixel(profile.indices[i].x(),profile.indices[i].y(),qRgb(0,255,0));
    }
    qDebug() << profile.indices.size() << profile.points.size() << "===============";
    imageWidget->setImage(curImage);
//    ui->label->setPixmap(pixmap->fromImage(curImage));
    qDebug() << 4;
    QVector<double> x, y;
    for(int i = 0; i < profile.points.size(); i++){
        x.push_back(profile.points[i].x());
        y.push_back(profile.points[i].y());
        qDebug() << profile.points[i];
    }
    GrWid->show();
    GrWid->plotGraph(x,y);
    GraphWidget* GrWid1 = new GraphWidget();
    GraphWidget* GrWid2 = new GraphWidget();
    QVector<QPointF> cutProf = ImageProcessor::cutProfile(profile.points, 5);
    QVector<double> x1, y1;
    QVector<double> x2, y2;
    for(int i = 0; i < cutProf.size(); i++){
        x1.push_back(cutProf[i].x());
        y1.push_back(cutProf[i].y());
    }
    GrWid1->plotGraph(x1,y1);
    GrWid1->show();
    edgeSelectionMode = false;
    qDebug() << 5;
    QVector<QPointF> toTestReInt = ImageProcessor::reInterpolateProfile(cutProf, 100);

    for(int i = 0; i < toTestReInt.size(); i++){
        x2.push_back(toTestReInt[i].x());
        y2.push_back(toTestReInt[i].y());
    }
    GrWid2->plotGraph(x2,y2);
    GrWid2->plotTwoGraphs(x1,y1,x2,y2);
    GrWid2->show();
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
//    setImage(sampleTwoHollows());
}

void MainWindow::on_actionTwo_hollows_big_triggered()
{
    QImage setImage = ImageProcessor::sampleTwoHollowsBig();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);
//    setImage(sampleTwoHollowsBig());
}

void MainWindow::on_actionimage_calculator_triggered()
{
    imCalculator->show();
}

void MainWindow::on_actionsharpen_triggered()
{
    Matrix2D<double> kernel = {{-1,-1,-1},
                               {-1, 9,-1},
                               {-1,-1,-1}};
    curImage = ImageProcessor::convImage(curImage, kernel);
//        ui->label->setPixmap(pixmap->fromImage(curImage));

    imageWidget->setImage(curImage);
    imageWidget->show();
    im1 = curImage;
}


