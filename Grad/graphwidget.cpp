// graphwidget.cpp
#include "graphwidget.h"
#include "ui_graphwidget.h"
#include <algorithm>

GraphWidget::GraphWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::GraphWidget)
{
    ui->setupUi(this);

    // Configure default plot settings / Nastraivayem parametry grafika po-umolchaniyu
    ui->widget->setBackground(Qt::white);
    ui->widget->axisRect()->setBackground(Qt::white);
    ui->widget->xAxis->setLabel("x");
    ui->widget->yAxis->setLabel("y");
    ui->widget->legend->setVisible(false);
}

GraphWidget::~GraphWidget()
{
    delete ui;
}

void GraphWidget::plotGraph(const QVector<double>& x, const QVector<double>& y)
{
    if (x.isEmpty() || y.isEmpty()) {
        qDebug() << "Error: Empty data in plotGraph";
        return;
    }

    if (x.size() != y.size()) {
        qDebug() << "Error: X and Y vectors have different sizes in plotGraph";
        return;
    }

    clearGraphs();

    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x, y);
    ui->widget->graph(0)->setPen(QPen(Qt::blue, 2));
    ui->widget->graph(0)->setName("Data");

    updateAxisRanges(x, y);
    ui->widget->replot();
}

void GraphWidget::plotTwoGraphs(const QVector<double>& x1, const QVector<double>& y1,
                                const QVector<double>& x2, const QVector<double>& y2)
{
    if (x1.isEmpty() || y1.isEmpty() || x2.isEmpty() || y2.isEmpty()) {
        qDebug() << "Error: Empty data in plotTwoGraphs";
        return;
    }

    if (x1.size() != y1.size()) {
        qDebug() << "Error: X1 and Y1 have different sizes";
        return;
    }

    if (x2.size() != y2.size()) {
        qDebug() << "Error: X2 and Y2 have different sizes";
        return;
    }

    clearGraphs();

    // First graph (blue) / Pervyy grafik (siniy)
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x1, y1);
    ui->widget->graph(0)->setPen(QPen(Qt::blue, 2));
    ui->widget->graph(0)->setName("Profile 1");

    // Second graph (red) / Vtoroy grafik (krasnyy)
    ui->widget->addGraph();
    ui->widget->graph(1)->setData(x2, y2);
    ui->widget->graph(1)->setPen(QPen(Qt::red, 2));
    ui->widget->graph(1)->setName("Profile 2");

    updateAxisRangesForTwoGraphs(x1, y1, x2, y2);

    ui->widget->legend->setVisible(true);
    ui->widget->replot();
}

void GraphWidget::clearGraphs()
{
    ui->widget->clearGraphs();
}

void GraphWidget::setAxisLabels(const QString& xLabel, const QString& yLabel)
{
    ui->widget->xAxis->setLabel(xLabel);
    ui->widget->yAxis->setLabel(yLabel);
}

void GraphWidget::setTitle(const QString& title)
{
    ui->widget->plotLayout()->insertRow(0);
    QCPTextElement* titleElement = new QCPTextElement(ui->widget, title, QFont("sans", 12, QFont::Bold));
    ui->widget->plotLayout()->addElement(0, 0, titleElement);
}

void GraphWidget::setLegendVisible(bool visible)
{
    ui->widget->legend->setVisible(visible);
}

void GraphWidget::updateAxisRanges(const QVector<double>& x, const QVector<double>& y)
{
    if (x.isEmpty() || y.isEmpty()) return;

    double minX = *std::min_element(x.begin(), x.end());
    double maxX = *std::max_element(x.begin(), x.end());
    double minY = *std::min_element(y.begin(), y.end());
    double maxY = *std::max_element(y.begin(), y.end());

    // Add 5% margin / Dobavlyayem 5% otstup
    double xRange = maxX - minX;
    double yRange = maxY - minY;

    if (xRange < 1e-10) xRange = 1.0;
    if (yRange < 1e-10) yRange = 1.0;

    ui->widget->xAxis->setRange(minX - 0.05 * xRange, maxX + 0.05 * xRange);
    ui->widget->yAxis->setRange(minY - 0.05 * yRange, maxY + 0.05 * yRange);
}

void GraphWidget::updateAxisRangesForTwoGraphs(const QVector<double>& x1, const QVector<double>& y1,
                                               const QVector<double>& x2, const QVector<double>& y2)
{
    if (x1.isEmpty() || y1.isEmpty() || x2.isEmpty() || y2.isEmpty()) return;

    double minX = std::min(*std::min_element(x1.begin(), x1.end()),
                           *std::min_element(x2.begin(), x2.end()));
    double maxX = std::max(*std::max_element(x1.begin(), x1.end()),
                           *std::max_element(x2.begin(), x2.end()));
    double minY = std::min(*std::min_element(y1.begin(), y1.end()),
                           *std::min_element(y2.begin(), y2.end()));
    double maxY = std::max(*std::max_element(y1.begin(), y1.end()),
                           *std::max_element(y2.begin(), y2.end()));

    // Add 5% margin / Dobavlyayem 5% otstup
    double xRange = maxX - minX;
    double yRange = maxY - minY;

    if (xRange < 1e-10) xRange = 1.0;
    if (yRange < 1e-10) yRange = 1.0;

    ui->widget->xAxis->setRange(minX - 0.05 * xRange, maxX + 0.05 * xRange);
    ui->widget->yAxis->setRange(minY - 0.05 * yRange, maxY + 0.05 * yRange);
}
