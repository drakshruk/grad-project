#include "imagecalculator.h"
#include "ui_imagecalculator.h"

ImageCalculator::ImageCalculator(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ImageCalculator)
{
    ui->setupUi(this);
    ui->image1ComboBox->addItem("Select other");
    ui->image2ComboBox->addItem("Select other");
    ui->operationComboBox->addItem("Add");
    ui->operationComboBox->addItem("Substract");
    ui->operationComboBox->addItem("Multiply");
    ui->operationComboBox->addItem("Divide");
    ui->operationComboBox->addItem("AND");
    ui->operationComboBox->addItem("OR");
    ui->operationComboBox->addItem("XOR");
    ui->operationComboBox->addItem("Min");
    ui->operationComboBox->addItem("Max");
    ui->operationComboBox->addItem("Average");
    ui->operationComboBox->addItem("Difference");
    ui->operationComboBox->addItem("Copy");
    ui->operationComboBox->addItem("Transparent-zero");
}

ImageCalculator::~ImageCalculator()
{
    delete ui;
}

void ImageCalculator::on_image1ComboBox_activated(int arg)
{
    Q_UNUSED(arg);
    if(ui->image1ComboBox->currentText() == "Select other"){
        QFileDialog* openImDialog = new QFileDialog(this);
        openImDialog->setFileMode(QFileDialog::AnyFile);
        openImDialog->setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));

        QString im1FileName = openImDialog->getOpenFileName();
        ui->image1ComboBox->addItem(im1FileName);
        ui->image1ComboBox->setCurrentIndex(ui->image1ComboBox->currentIndex() + 1);
        image1 = QImage(im1FileName);
    }
}


void ImageCalculator::on_image2ComboBox_activated(int arg)
{
    Q_UNUSED(arg);
    if(ui->image2ComboBox->currentText() == "Select other"){
        QFileDialog* openImDialog = new QFileDialog(this);
        openImDialog->setFileMode(QFileDialog::AnyFile);
        openImDialog->setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));

        QString im2FileName = openImDialog->getOpenFileName();
        ui->image2ComboBox->addItem(im2FileName);
        ui->image2ComboBox->setCurrentIndex(ui->image2ComboBox->currentIndex() + 1);
        image2 = QImage(im2FileName);
    }
}


void ImageCalculator::on_operationComboBox_activated(int arg)
{
    Q_UNUSED(arg);
    operation = ui->operationComboBox->currentText();
}


void ImageCalculator::on_newWindowRadioButton_clicked(bool checked)
{
    newWindow = checked;
}


void ImageCalculator::on_floatResRadioButton_clicked(bool checked)
{
    floatResult = checked;
}


void ImageCalculator::on_calculateButton_clicked()
{
    emit throw_imageCalculator(image1, image2, operation, newWindow, floatResult);
    this->close();
}


void ImageCalculator::on_closeButton_clicked()
{
    this->close();
}

