/********************************************************************************
** Form generated from reading UI file 'graphwidget.ui'
**
** Created by: Qt User Interface Compiler version 5.12.12
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GRAPHWIDGET_H
#define UI_GRAPHWIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_GraphWidget
{
public:
    QGridLayout *gridLayout;
    QCustomPlot *widget;

    void setupUi(QWidget *GraphWidget)
    {
        if (GraphWidget->objectName().isEmpty())
            GraphWidget->setObjectName(QString::fromUtf8("GraphWidget"));
        GraphWidget->resize(400, 300);
        gridLayout = new QGridLayout(GraphWidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        widget = new QCustomPlot(GraphWidget);
        widget->setObjectName(QString::fromUtf8("widget"));

        gridLayout->addWidget(widget, 0, 0, 1, 1);


        retranslateUi(GraphWidget);

        QMetaObject::connectSlotsByName(GraphWidget);
    } // setupUi

    void retranslateUi(QWidget *GraphWidget)
    {
        GraphWidget->setWindowTitle(QApplication::translate("GraphWidget", "Form", nullptr));
    } // retranslateUi

};

namespace Ui {
    class GraphWidget: public Ui_GraphWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GRAPHWIDGET_H
