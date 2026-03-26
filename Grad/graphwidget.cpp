#include "graphwidget.h"
#include "ui_graphwidget.h"

GraphWidget::GraphWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::GraphWidget)
{
    ui->setupUi(this);
}

GraphWidget::~GraphWidget()
{
    delete ui;
}

void GraphWidget::plotGraph(QVector<double> x, QVector<double> y)
{
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x, y);
    ui->widget->xAxis->setLabel("x");
    ui->widget->yAxis->setLabel("y");
    ui->widget->xAxis->setRange(x[0],x.back());
    double minY = y[0];
    for(int i = 0; i < y.size(); i++)
    {
        if(minY > y[i]) minY = y[i];
    }
    double maxY = y[0];
    for(int i = 0; i < y.size(); i++)
    {
        if(maxY < y[i]) maxY = y[i];
    }
    ui->widget->yAxis->setRange(minY,maxY);
    ui->widget->replot();
}
void GraphWidget::plotTwoGraphs(QVector<double> x1, QVector<double> y1, QVector<double> x2, QVector<double> y2)
{
    ui->widget->clearGraphs();

    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x1, y1);
    ui->widget->graph(0)->setPen(QPen(Qt::blue, 2)); // Синий цвет, толщина 2
    ui->widget->graph(0)->setName("Profile 1"); // Легенда

    ui->widget->addGraph();
    ui->widget->graph(1)->setData(x2, y2);
    ui->widget->graph(1)->setPen(QPen(Qt::red, 2)); // Красный цвет, толщина 2
    ui->widget->graph(1)->setName("Profile 2"); // Легенда

    // Настройка осей
    ui->widget->xAxis->setLabel("Position");
    ui->widget->yAxis->setLabel("Intensity");

    if (x1.isEmpty() || y1.isEmpty() || x2.isEmpty() || y2.isEmpty()) {
        qDebug() << "Error: Empty data vectors";
        return;
    }

    double minX = std::min(*std::min_element(x1.begin(), x1.end()),
                          *std::min_element(x2.begin(), x2.end()));
    double maxX = std::max(*std::max_element(x1.begin(), x1.end()),
                          *std::max_element(x2.begin(), x2.end()));

    double minY = std::min(*std::min_element(y1.begin(), y1.end()),
                          *std::min_element(y2.begin(), y2.end()));
    double maxY = std::max(*std::max_element(y1.begin(), y1.end()),
                          *std::max_element(y2.begin(), y2.end()));

    double xRange = maxX - minX;
    double yRange = maxY - minY;

    ui->widget->xAxis->setRange(minX - 0.05 * xRange, maxX + 0.05 * xRange);
    ui->widget->yAxis->setRange(minY - 0.05 * yRange, maxY + 0.05 * yRange);

    ui->widget->legend->setVisible(true);

    ui->widget->replot();

//    qDebug() << "Plotted two graphs:";
//    qDebug() << "Graph 1 - Points:" << x1.size() << "X range: [" << *std::min_element(x1.begin(), x1.end())
//             << "," << *std::max_element(x1.begin(), x1.end()) << "]";
//    qDebug() << "Graph 2 - Points:" << x2.size() << "X range: [" << *std::min_element(x2.begin(), x2.end())
//             << "," << *std::max_element(x2.begin(), x2.end()) << "]";
}
