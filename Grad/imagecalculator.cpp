// imagecalculator.cpp
#include "imagecalculator.h"
#include "ui_imagecalculator.h"

// Helper function to convert operation to string / Vspomogatel'naya funktsiya dlya preobrazovaniya operatsii v stroku
QString ImageCalculator::operationToString(ImageOperation op)
{
    switch (op) {
    case ImageOperation::Add: return "Add";
    case ImageOperation::Subtract: return "Subtract";
    case ImageOperation::Multiply: return "Multiply";
    case ImageOperation::Divide: return "Divide";
    case ImageOperation::AND: return "AND";
    case ImageOperation::OR: return "OR";
    case ImageOperation::XOR: return "XOR";
    case ImageOperation::Min: return "Min";
    case ImageOperation::Max: return "Max";
    case ImageOperation::Average: return "Average";
    case ImageOperation::Difference: return "Difference";
    case ImageOperation::Copy: return "Copy";
    case ImageOperation::TransparentZero: return "Transparent-zero";
    default: return "Add";
    }
}

// Helper function to convert string to operation / Vspomogatel'naya funktsiya dlya preobrazovaniya stroki v operatsiyu
ImageOperation ImageCalculator::stringToOperation(const QString& str)
{
    if (str == "Add") return ImageOperation::Add;
    if (str == "Subtract") return ImageOperation::Subtract;
    if (str == "Multiply") return ImageOperation::Multiply;
    if (str == "Divide") return ImageOperation::Divide;
    if (str == "AND") return ImageOperation::AND;
    if (str == "OR") return ImageOperation::OR;
    if (str == "XOR") return ImageOperation::XOR;
    if (str == "Min") return ImageOperation::Min;
    if (str == "Max") return ImageOperation::Max;
    if (str == "Average") return ImageOperation::Average;
    if (str == "Difference") return ImageOperation::Difference;
    if (str == "Copy") return ImageOperation::Copy;
    if (str == "Transparent-zero") return ImageOperation::TransparentZero;
    return ImageOperation::Add;
}

ImageCalculator::ImageCalculator(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ImageCalculator)
{
    ui->setupUi(this);

    // Populate combo boxes / Zapolnyayem combo boxes
    ui->image1ComboBox->addItem("Select other");
    ui->image2ComboBox->addItem("Select other");

    populateOperationComboBox();

    // Set default selection / Ustanavlivaem vybor po-umolchaniyu
    ui->newWindowRadioButton->setChecked(true);
}

ImageCalculator::~ImageCalculator()
{
    delete ui;
}

void ImageCalculator::populateOperationComboBox()
{
    // Add all operations to combo box / Dobavlyayem vse operatsii v combo box
    ui->operationComboBox->addItem("Add");
    ui->operationComboBox->addItem("Subtract");
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

    // Set default / Ustanavlivaem po-umolchaniyu
    ui->operationComboBox->setCurrentIndex(0);
    m_currentOperation = ImageOperation::Add;
}

void ImageCalculator::setImage1(const QImage& image, const QString& name)
{
    m_image1 = image;
    if (!name.isEmpty()) {
        ui->image1ComboBox->addItem(name);
        ui->image1ComboBox->setCurrentIndex(ui->image1ComboBox->count() - 1);
    }
}

void ImageCalculator::setImage2(const QImage& image, const QString& name)
{
    m_image2 = image;
    if (!name.isEmpty()) {
        ui->image2ComboBox->addItem(name);
        ui->image2ComboBox->setCurrentIndex(ui->image2ComboBox->count() - 1);
    }
}

void ImageCalculator::loadImageFromFile(int imageIndex)
{
    QFileDialog openImDialog(this);
    openImDialog.setFileMode(QFileDialog::AnyFile);
    openImDialog.setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));

    QString fileName = openImDialog.getOpenFileName();
    if (fileName.isEmpty()) return;

    QImage loadedImage(fileName);
    if (loadedImage.isNull()) {
        qDebug() << "Failed to load image:" << fileName;
        return;
    }

    // Store the image and update combo box / Sohranyayem izobrazheniye i obnovlyayem combo box
    if (imageIndex == 1) {
        m_image1 = loadedImage;
        ui->image1ComboBox->addItem(fileName);
        ui->image1ComboBox->setCurrentIndex(ui->image1ComboBox->count() - 1);
    } else {
        m_image2 = loadedImage;
        ui->image2ComboBox->addItem(fileName);
        ui->image2ComboBox->setCurrentIndex(ui->image2ComboBox->count() - 1);
    }
}

void ImageCalculator::on_image1ComboBox_activated(int index)
{
    Q_UNUSED(index);
    if (ui->image1ComboBox->currentText() == "Select other") {
        loadImageFromFile(1);
    }
}

void ImageCalculator::on_image2ComboBox_activated(int index)
{
    Q_UNUSED(index);
    if (ui->image2ComboBox->currentText() == "Select other") {
        loadImageFromFile(2);
    }
}

void ImageCalculator::on_operationComboBox_activated(int index)
{
    Q_UNUSED(index);
    m_currentOperation = stringToOperation(ui->operationComboBox->currentText());
}

void ImageCalculator::on_newWindowRadioButton_clicked(bool checked)
{
    m_newWindow = checked;
}

void ImageCalculator::on_floatResRadioButton_clicked(bool checked)
{
    m_floatResult = checked;
}

void ImageCalculator::on_calculateButton_clicked()
{
    // Validate images / Proveryayem izobrazheniya
    if (m_image1.isNull()) {
        qDebug() << "Error: Image 1 is not loaded";
        return;
    }

    if (m_image2.isNull() && m_currentOperation != ImageOperation::Copy) {
        qDebug() << "Error: Image 2 is not loaded";
        return;
    }

    // Emit both signals for compatibility / Ispuskayem oba signala dlya sovmestimosti
    QString opString = operationToString(m_currentOperation);
    emit calculateRequested(m_image1, m_image2, m_currentOperation, m_newWindow, m_floatResult);
    emit throw_imageCalculator(m_image1, m_image2, opString, m_newWindow, m_floatResult);

    this->close();
}

void ImageCalculator::on_closeButton_clicked()
{
    this->close();
}
