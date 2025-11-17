/********************************************************************************
** Form generated from reading UI file 'imageshowcasewidget.ui'
**
** Created by: Qt User Interface Compiler version 5.12.12
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_IMAGESHOWCASEWIDGET_H
#define UI_IMAGESHOWCASEWIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ImageShowcaseWidget
{
public:

    void setupUi(QWidget *ImageShowcaseWidget)
    {
        if (ImageShowcaseWidget->objectName().isEmpty())
            ImageShowcaseWidget->setObjectName(QString::fromUtf8("ImageShowcaseWidget"));
        ImageShowcaseWidget->resize(400, 300);

        retranslateUi(ImageShowcaseWidget);

        QMetaObject::connectSlotsByName(ImageShowcaseWidget);
    } // setupUi

    void retranslateUi(QWidget *ImageShowcaseWidget)
    {
        ImageShowcaseWidget->setWindowTitle(QApplication::translate("ImageShowcaseWidget", "Form", nullptr));
    } // retranslateUi

};

namespace Ui {
    class ImageShowcaseWidget: public Ui_ImageShowcaseWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IMAGESHOWCASEWIDGET_H
