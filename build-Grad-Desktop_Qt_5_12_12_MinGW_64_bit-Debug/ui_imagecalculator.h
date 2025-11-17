/********************************************************************************
** Form generated from reading UI file 'imagecalculator.ui'
**
** Created by: Qt User Interface Compiler version 5.12.12
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_IMAGECALCULATOR_H
#define UI_IMAGECALCULATOR_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ImageCalculator
{
public:
    QGridLayout *gridLayout;
    QComboBox *image1ComboBox;
    QPushButton *closeButton;
    QLabel *label_3;
    QComboBox *operationComboBox;
    QLabel *label;
    QComboBox *image2ComboBox;
    QRadioButton *floatResRadioButton;
    QPushButton *calculateButton;
    QLabel *label_2;
    QRadioButton *newWindowRadioButton;

    void setupUi(QWidget *ImageCalculator)
    {
        if (ImageCalculator->objectName().isEmpty())
            ImageCalculator->setObjectName(QString::fromUtf8("ImageCalculator"));
        ImageCalculator->resize(412, 458);
        gridLayout = new QGridLayout(ImageCalculator);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        image1ComboBox = new QComboBox(ImageCalculator);
        image1ComboBox->setObjectName(QString::fromUtf8("image1ComboBox"));

        gridLayout->addWidget(image1ComboBox, 0, 1, 1, 3);

        closeButton = new QPushButton(ImageCalculator);
        closeButton->setObjectName(QString::fromUtf8("closeButton"));

        gridLayout->addWidget(closeButton, 5, 2, 1, 1);

        label_3 = new QLabel(ImageCalculator);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setMaximumSize(QSize(75, 24));

        gridLayout->addWidget(label_3, 2, 0, 1, 1);

        operationComboBox = new QComboBox(ImageCalculator);
        operationComboBox->setObjectName(QString::fromUtf8("operationComboBox"));

        gridLayout->addWidget(operationComboBox, 1, 1, 1, 3);

        label = new QLabel(ImageCalculator);
        label->setObjectName(QString::fromUtf8("label"));
        label->setMaximumSize(QSize(75, 24));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        image2ComboBox = new QComboBox(ImageCalculator);
        image2ComboBox->setObjectName(QString::fromUtf8("image2ComboBox"));

        gridLayout->addWidget(image2ComboBox, 2, 1, 1, 3);

        floatResRadioButton = new QRadioButton(ImageCalculator);
        floatResRadioButton->setObjectName(QString::fromUtf8("floatResRadioButton"));

        gridLayout->addWidget(floatResRadioButton, 4, 0, 1, 3);

        calculateButton = new QPushButton(ImageCalculator);
        calculateButton->setObjectName(QString::fromUtf8("calculateButton"));

        gridLayout->addWidget(calculateButton, 5, 1, 1, 1);

        label_2 = new QLabel(ImageCalculator);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setMaximumSize(QSize(75, 24));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        newWindowRadioButton = new QRadioButton(ImageCalculator);
        newWindowRadioButton->setObjectName(QString::fromUtf8("newWindowRadioButton"));
        newWindowRadioButton->setChecked(true);

        gridLayout->addWidget(newWindowRadioButton, 3, 0, 1, 3);


        retranslateUi(ImageCalculator);

        QMetaObject::connectSlotsByName(ImageCalculator);
    } // setupUi

    void retranslateUi(QWidget *ImageCalculator)
    {
        ImageCalculator->setWindowTitle(QApplication::translate("ImageCalculator", "Form", nullptr));
        closeButton->setText(QApplication::translate("ImageCalculator", "Close", nullptr));
        label_3->setText(QApplication::translate("ImageCalculator", "Image 2:", nullptr));
        label->setText(QApplication::translate("ImageCalculator", "Image 1:", nullptr));
        floatResRadioButton->setText(QApplication::translate("ImageCalculator", "32-bit (float) result", nullptr));
        calculateButton->setText(QApplication::translate("ImageCalculator", "Ok", nullptr));
        label_2->setText(QApplication::translate("ImageCalculator", "Operation:", nullptr));
        newWindowRadioButton->setText(QApplication::translate("ImageCalculator", "Create new window", nullptr));
    } // retranslateUi

};

namespace Ui {
    class ImageCalculator: public Ui_ImageCalculator {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IMAGECALCULATOR_H
