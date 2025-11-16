#ifndef GRAPHWIDGET_H
#define GRAPHWIDGET_H

#include <QWidget>
#include "qcustomplot.h"

namespace Ui {
class GraphWidget;
}

class GraphWidget : public QWidget
{
    Q_OBJECT

public:
    explicit GraphWidget(QWidget *parent = nullptr);
    ~GraphWidget();
    void plotGraph(QVector<double> x, QVector<double> y);
    void plotTwoGraphs(QVector<double> x1, QVector<double> y1, QVector<double> x2, QVector<double> y2);

private:
    Ui::GraphWidget *ui;
};

#endif // GRAPHWIDGET_H
