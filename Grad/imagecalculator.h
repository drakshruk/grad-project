// imagecalculator.h
#ifndef IMAGECALCULATOR_H
#define IMAGECALCULATOR_H

#include <QWidget>
#include <QImage>
#include <QFileDialog>
#include <QDebug>

namespace Ui {
class ImageCalculator;
}

/**
 * EN: Image operation types for calculator
 * RU: Tipy operatsiy dlya kalkulyatora izobrazheniy
 */
enum class ImageOperation {
    Add,                // Addition / Slozheniye
    Subtract,           // Subtraction / Vychitaniye
    Multiply,           // Multiplication / Umnozheniye
    Divide,             // Division / Deleniye
    AND,                // Bitwise AND / Pobitovoye I
    OR,                 // Bitwise OR / Pobitovoye ILI
    XOR,                // Bitwise XOR / Pobitovoye isklyuchayushcheye ILI
    Min,                // Minimum / Minimum
    Max,                // Maximum / Maximum
    Average,            // Average / Sredneye
    Difference,         // Absolute difference / Modul' raznosti
    Copy,               // Copy first image / Kopirovat' pervoye izobrazheniye
    TransparentZero     // Make zero pixels transparent / Sdelat' nulevyye pikseli prozrachnymi
};

/**
 * EN: Image calculator widget for performing pixel-wise operations between two images
 * RU: Vidzhet kalkulyatora izobrazheniy dlya vypolneniya popelementnykh operatsiy mezhdu dvumya izobrazheniyami
 */
class ImageCalculator : public QWidget
{
    Q_OBJECT

public:
    explicit ImageCalculator(QWidget *parent = nullptr);
    ~ImageCalculator();

    /**
     * EN: Sets the first image (from file or current)
     * RU: Ustanavlivayet pervoye izobrazheniye (iz fayla ili tekushcheye)
     */
    void setImage1(const QImage& image, const QString& name = "");

    /**
     * EN: Sets the second image
     * RU: Ustanavlivayet vtoroye izobrazheniye
     */
    void setImage2(const QImage& image, const QString& name = "");

    /**
     * EN: Gets the currently selected operation
     * RU: Poluchayet vybrannuyu operatsiyu
     */
    ImageOperation getCurrentOperation() const { return m_currentOperation; }

signals:
    /**
     * EN: Emitted when calculation is requested with all parameters
     * RU: Ispuskayetsya pri zaprose vychisleniya so vsemi parametrami
     */
    void calculateRequested(const QImage& image1, const QImage& image2,
                            ImageOperation operation, bool newWindow, bool floatResult);

    /**
     * EN: Backward compatibility signal for existing code
     * RU: Signal obratnoy sovmestimosti dlya sushchestvuyushchego koda
     */
    void throw_imageCalculator(const QImage& image1, const QImage& image2,
                               QString operation, bool newWindow, bool floatResult);

private slots:
    void on_image1ComboBox_activated(int index);
    void on_image2ComboBox_activated(int index);
    void on_operationComboBox_activated(int index);
    void on_newWindowRadioButton_clicked(bool checked);
    void on_floatResRadioButton_clicked(bool checked);
    void on_calculateButton_clicked();
    void on_closeButton_clicked();

private:
    /**
     * EN: Converts operation enum to string for backward compatibility
     * RU: Preobrazuyet enum operatsii v stroku dlya obratnoy sovmestimosti
     */
    QString operationToString(ImageOperation op);

    /**
     * EN: Converts string to operation enum
     * RU: Preobrazuyet stroku v enum operatsii
     */
    ImageOperation stringToOperation(const QString& str);

    /**
     * EN: Populates operation combo box with all operations
     * RU: Zapolnyayet combo box operatsiyami
     */
    void populateOperationComboBox();

    /**
     * EN: Loads image from file and adds to combo box
     * RU: Zagruzhaet izobrazheniye iz fayla i dobavlyayet v combo box
     */
    void loadImageFromFile(int imageIndex); // 1 or 2

private:
    QImage m_image1;
    QImage m_image2;
    ImageOperation m_currentOperation = ImageOperation::Add;
    bool m_newWindow = true;
    bool m_floatResult = false;
    Ui::ImageCalculator *ui;
};

#endif // IMAGECALCULATOR_H
