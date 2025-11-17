/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.12.12
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actiongaussian_blur;
    QAction *actiongaussian_edge_detection;
    QAction *actionpcr_edge_detection;
    QAction *actionopen_file;
    QAction *actionsave_file;
    QAction *actiongradient_X_and_Y;
    QAction *actiontest;
    QAction *actiondraw_profile;
    QAction *actionTwo_hollows;
    QAction *actionTwo_hollows_big;
    QAction *actionLaplacian_edge_detection;
    QAction *actionimage_calculator;
    QAction *actionsharpen;
    QWidget *centralWidget;
    QGridLayout *gridLayout;
    QLabel *label;
    QMenuBar *menuBar;
    QMenu *menufile;
    QMenu *menumath;
    QMenu *menutools;
    QMenu *menusamples;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(658, 107);
        actiongaussian_blur = new QAction(MainWindow);
        actiongaussian_blur->setObjectName(QString::fromUtf8("actiongaussian_blur"));
        actiongaussian_edge_detection = new QAction(MainWindow);
        actiongaussian_edge_detection->setObjectName(QString::fromUtf8("actiongaussian_edge_detection"));
        actionpcr_edge_detection = new QAction(MainWindow);
        actionpcr_edge_detection->setObjectName(QString::fromUtf8("actionpcr_edge_detection"));
        actionopen_file = new QAction(MainWindow);
        actionopen_file->setObjectName(QString::fromUtf8("actionopen_file"));
        actionsave_file = new QAction(MainWindow);
        actionsave_file->setObjectName(QString::fromUtf8("actionsave_file"));
        actiongradient_X_and_Y = new QAction(MainWindow);
        actiongradient_X_and_Y->setObjectName(QString::fromUtf8("actiongradient_X_and_Y"));
        actiontest = new QAction(MainWindow);
        actiontest->setObjectName(QString::fromUtf8("actiontest"));
        actiondraw_profile = new QAction(MainWindow);
        actiondraw_profile->setObjectName(QString::fromUtf8("actiondraw_profile"));
        actionTwo_hollows = new QAction(MainWindow);
        actionTwo_hollows->setObjectName(QString::fromUtf8("actionTwo_hollows"));
        actionTwo_hollows_big = new QAction(MainWindow);
        actionTwo_hollows_big->setObjectName(QString::fromUtf8("actionTwo_hollows_big"));
        actionLaplacian_edge_detection = new QAction(MainWindow);
        actionLaplacian_edge_detection->setObjectName(QString::fromUtf8("actionLaplacian_edge_detection"));
        actionimage_calculator = new QAction(MainWindow);
        actionimage_calculator->setObjectName(QString::fromUtf8("actionimage_calculator"));
        actionsharpen = new QAction(MainWindow);
        actionsharpen->setObjectName(QString::fromUtf8("actionsharpen"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        gridLayout = new QGridLayout(centralWidget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label = new QLabel(centralWidget);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 658, 25));
        menufile = new QMenu(menuBar);
        menufile->setObjectName(QString::fromUtf8("menufile"));
        menumath = new QMenu(menuBar);
        menumath->setObjectName(QString::fromUtf8("menumath"));
        menutools = new QMenu(menuBar);
        menutools->setObjectName(QString::fromUtf8("menutools"));
        menusamples = new QMenu(menuBar);
        menusamples->setObjectName(QString::fromUtf8("menusamples"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menufile->menuAction());
        menuBar->addAction(menumath->menuAction());
        menuBar->addAction(menutools->menuAction());
        menuBar->addAction(menusamples->menuAction());
        menufile->addAction(actionopen_file);
        menufile->addAction(actionsave_file);
        menumath->addAction(actionsharpen);
        menumath->addAction(actiongaussian_blur);
        menumath->addAction(actiongaussian_edge_detection);
        menumath->addAction(actionpcr_edge_detection);
        menumath->addAction(actiongradient_X_and_Y);
        menumath->addAction(actionLaplacian_edge_detection);
        menutools->addAction(actiontest);
        menutools->addAction(actiondraw_profile);
        menutools->addAction(actionimage_calculator);
        menusamples->addAction(actionTwo_hollows);
        menusamples->addAction(actionTwo_hollows_big);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", nullptr));
        actiongaussian_blur->setText(QApplication::translate("MainWindow", "gaussian blur", nullptr));
        actiongaussian_edge_detection->setText(QApplication::translate("MainWindow", "gaussian edge detection", nullptr));
        actionpcr_edge_detection->setText(QApplication::translate("MainWindow", "pcr edge detection", nullptr));
        actionopen_file->setText(QApplication::translate("MainWindow", "open file", nullptr));
        actionsave_file->setText(QApplication::translate("MainWindow", "save file", nullptr));
        actiongradient_X_and_Y->setText(QApplication::translate("MainWindow", "gradient X and Y", nullptr));
        actiontest->setText(QApplication::translate("MainWindow", "test", nullptr));
        actiondraw_profile->setText(QApplication::translate("MainWindow", "draw profile", nullptr));
        actionTwo_hollows->setText(QApplication::translate("MainWindow", "two hollows", nullptr));
        actionTwo_hollows_big->setText(QApplication::translate("MainWindow", "two hollows big", nullptr));
        actionLaplacian_edge_detection->setText(QApplication::translate("MainWindow", "Laplacian edge detection", nullptr));
        actionimage_calculator->setText(QApplication::translate("MainWindow", "image calculator", nullptr));
        actionsharpen->setText(QApplication::translate("MainWindow", "sharpen", nullptr));
        label->setText(QString());
        menufile->setTitle(QApplication::translate("MainWindow", "file", nullptr));
        menumath->setTitle(QApplication::translate("MainWindow", "math", nullptr));
        menutools->setTitle(QApplication::translate("MainWindow", "tools", nullptr));
        menusamples->setTitle(QApplication::translate("MainWindow", "samples", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
