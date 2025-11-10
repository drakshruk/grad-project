#ifndef IMAGECALCULATOR_H
#define IMAGECALCULATOR_H

#include <QWidget>
#include <QFileDialog>

namespace Ui {
class ImageCalculator;
}

class ImageCalculator : public QWidget
{
    Q_OBJECT

public:
    explicit ImageCalculator(QWidget *parent = nullptr);
    ~ImageCalculator();

private slots:
    void on_image1ComboBox_activated(int arg);

    void on_image2ComboBox_activated(int arg);

    void on_operationComboBox_activated(int arg);

    void on_newWindowRadioButton_clicked(bool checked);

    void on_floatResRadioButton_clicked(bool checked);

    void on_calculateButton_clicked();

    void on_closeButton_clicked();

signals:
    void throw_imageCalculator(const QImage& image1, const QImage& image2, QString operation, bool newWindow, bool floatResult);

private:
    QImage image1, image2;
    QString operation;
    bool newWindow = true, floatResult = false;
    Ui::ImageCalculator *ui;
};

#endif // IMAGECALCULATOR_H
