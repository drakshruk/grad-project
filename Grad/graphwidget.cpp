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
    ui->widget->addGraph();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x1, y1);
    ui->widget->graph(1)->setData(x2, y2);
    ui->widget->xAxis->setLabel("x");
    ui->widget->yAxis->setLabel("y");
    ui->widget->xAxis->setRange(x1[0],x1.back());
    double minY = y1[0];
    double minX = x1[0];
    for(int i = 0; i < y1.size(); i++)
    {
        if(minY > y1[i]) minY = y1[i];
        if(minX > x1[i]) minX = x1[i];

    }
    for(int i = 0; i < y2.size(); i++)
    {
        if(minY > y2[i]) minY = y2[i];
        if(minX > x2[i]) minX = x2[i];
    }
    double maxY = y1[0];
    double maxX = x1[0];
    for(int i = 0; i < y1.size(); i++)
    {
        if(maxY < y1[i]) maxY = y1[i];
        if(maxX < x1[i]) maxX = x1[i];
    }
    for(int i = 0; i < y2.size(); i++)
    {
        if(maxY < y2[i]) maxY = y2[i];
        if(maxX < x2[i]) maxX = x2[i];
    }
    ui->widget->yAxis->setRange(minY,maxY);
    ui->widget->xAxis->setRange(minX,maxX);
    ui->widget->replot();
}
