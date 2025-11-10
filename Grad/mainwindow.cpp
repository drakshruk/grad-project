#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "utility.cpp"

/*
 *
 */
QPoint mPos;

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

    iRad = 10;
    dSigma = 5;
    dRedMax = 255;
    dBlueMax = 255;

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
        dSigma = getDataDialog->value;
        iRad = 3*dSigma;
        curImage = ImageProcessor::convImage(curImage, ImageProcessor::getGauss(iRad,iRad,dSigma));
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
        QVector<QPoint> selectedEdge;
        ImageProcessor::selectEdge(mPos, attMat, selectedEdge);
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
    ProfileResult cutProf = ImageProcessor::cutProfile(profile, 5);
    QVector<double> x1, y1;
    QVector<double> x2, y2;
    for(int i = 0; i < cutProf.points.size(); i++){
        x1.push_back(cutProf.points[i].x());
        y1.push_back(cutProf.points[i].y());
    }
    GrWid1->plotGraph(x1,y1);
    GrWid1->show();
    edgeSelectionMode = false;
    qDebug() << 5;
    ProfileResult toTestReInt = ImageProcessor::reInterpolateProfile(cutProf, 100);

    for(int i = 0; i < toTestReInt.points.size(); i++){
        x2.push_back(toTestReInt.points[i].x());
        y2.push_back(toTestReInt.points[i].y());
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
    m_dataModel->setCurrentImage(ImageProcessor::sampleTwoHollows());
//    setImage(sampleTwoHollows());
}

void MainWindow::on_actionTwo_hollows_big_triggered()
{
    m_dataModel->setCurrentImage(ImageProcessor::sampleTwoHollowsBig());
//    setImage(sampleTwoHollowsBig());
}

void MainWindow::on_actionimage_calculator_triggered()
{
    imCalculator->show();
}

void MainWindow::on_actionsharpen_triggered()
{
    Matrix2D<double> kernel = {{-1,-1,-1},
                               {-1,12,-1},
                               {-1,-1,-1}};
    curImage = ImageProcessor::convImage(curImage, kernel);
//        ui->label->setPixmap(pixmap->fromImage(curImage));

    imageWidget->setImage(curImage);
    imageWidget->show();
    im1 = curImage;
}


