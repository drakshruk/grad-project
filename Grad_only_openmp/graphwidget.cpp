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


    // ========================================================================
    // НАСТРОЙКА МАСШТАБИРОВАНИЯ ШРИФТОВ
    // ========================================================================

    // Устанавливаем размеры шрифтов
    QFont axisFont = ui->widget->xAxis->labelFont();
    axisFont.setPointSize(10);
    ui->widget->xAxis->setLabelFont(axisFont);
    ui->widget->yAxis->setLabelFont(axisFont);

    QFont tickFont = ui->widget->xAxis->tickLabelFont();
    tickFont.setPointSize(9);
    ui->widget->xAxis->setTickLabelFont(tickFont);
    ui->widget->yAxis->setTickLabelFont(tickFont);

    QFont legendFont = ui->widget->legend->font();
    legendFont.setPointSize(9);
    ui->widget->legend->setFont(legendFont);

    // Включаем автоматическое масштабирование шрифтов при изменении размера
    ui->widget->setAutoAddPlottableToLegend(true);
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

void GraphWidget::plotMultipleProfiles(const QVector<QVector<QPointF>>& profiles,
                                       const QVector<QString>& labels)
{
    if (profiles.isEmpty()) {
        qDebug() << "Error: Empty profiles data";
        return;
    }

    clearGraphs();

    // Цвета для разных кривых
    QVector<QColor> colors = {Qt::black, Qt::blue, Qt::green, Qt::red};

    for (int i = 0; i < profiles.size(); i++) {
        if (profiles[i].isEmpty()) continue;

        // Извлекаем X и Y координаты
        QVector<double> x(profiles[i].size());
        QVector<double> y(profiles[i].size());
        for (int j = 0; j < profiles[i].size(); j++) {
            x[j] = profiles[i][j].x();
            y[j] = profiles[i][j].y();
        }

        ui->widget->addGraph();
        ui->widget->graph(i)->setData(x, y);
        ui->widget->graph(i)->setPen(QPen(colors[i % colors.size()], 2));

        if (i < labels.size()) {
            ui->widget->graph(i)->setName(labels[i]);
        } else {
            ui->widget->graph(i)->setName(QString("Profile %1").arg(i));
        }
    }

    // Настройка осей
    ui->widget->xAxis->setLabel("X coordinate (pixels)");
    ui->widget->yAxis->setLabel("Intensity (0-255)");

    // Автоматический диапазон с отступом 5%
    ui->widget->rescaleAxes(true);
    ui->widget->xAxis->setRange(0, 200);
    ui->widget->yAxis->setRange(-10, 265);

    ui->widget->legend->setVisible(true);
    ui->widget->replot();
}

void GraphWidget::plotMultipleProfilesWithZeroCrossings(const QVector<QVector<QPointF>>& profiles,
                                                        const QVector<QString>& labels,
                                                        const QVector<QVector<double>>& zeroCrossings)
{
    if (profiles.isEmpty()) {
        qDebug() << "Error: Empty profiles data";
        return;
    }

    clearGraphs();

    // Цвета для разных кривых
    QVector<QColor> colors = {Qt::black, Qt::blue, Qt::green, Qt::red, Qt::magenta, Qt::cyan};

    for (int i = 0; i < profiles.size(); i++) {
        if (profiles[i].isEmpty()) continue;

        // Извлекаем X и Y координаты
        QVector<double> x(profiles[i].size());
        QVector<double> y(profiles[i].size());
        for (int j = 0; j < profiles[i].size(); j++) {
            x[j] = profiles[i][j].x();
            y[j] = profiles[i][j].y();
        }

        ui->widget->addGraph();
        ui->widget->graph(i)->setData(x, y);
        ui->widget->graph(i)->setPen(QPen(colors[i % colors.size()], 2));

        if (i < labels.size()) {
            ui->widget->graph(i)->setName(labels[i]);
        } else {
            ui->widget->graph(i)->setName(QString("Profile %1").arg(i));
        }

        // Отмечаем точки пересечения с нулём
        if (i < zeroCrossings.size()) {
            for (double zeroX : zeroCrossings[i]) {
                // Находим соответствующее значение Y (должно быть близко к 0)
                double zeroY = 0.0;
                // Добавляем вертикальную линию в точке zeroX
                QCPItemStraightLine* line = new QCPItemStraightLine(ui->widget);
                line->point1->setCoords(zeroX, -1000);
                line->point2->setCoords(zeroX, 1000);
                line->setPen(QPen(Qt::red, 1, Qt::DashLine));
            }
        }
    }

    // Настройка осей
    ui->widget->xAxis->setLabel("Distance along profile (pixels)");
    ui->widget->yAxis->setLabel("Laplacian response");

    // Автоматический диапазон с отступом 10%
    ui->widget->rescaleAxes(true);

    ui->widget->legend->setVisible(true);
    ui->widget->replot();
}

QCustomPlot* GraphWidget::getPlot()
{
    return ui->widget;
}