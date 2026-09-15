#ifndef GRAPHWIDGET_H
#define GRAPHWIDGET_H

#include <QWidget>
#include <QVector>
#include <QDebug>
#include "qcustomplot.h"

namespace Ui {
class GraphWidget;
}

/**
 * EN: Widget for displaying graphs (profiles, etc.)
 * RU: Vidzhet dlya otobrazheniya grafikov (profiley, i t.d.)
 */
class GraphWidget : public QWidget
{
    Q_OBJECT

public:
    explicit GraphWidget(QWidget *parent = nullptr);
    ~GraphWidget();

    /**
     * EN: Plots a single graph with given x and y data
     * RU: Otobrazhayet odin grafik s zadannymi dannymi x i y
     */
    void plotGraph(const QVector<double>& x, const QVector<double>& y);

    /**
     * EN: Plots two graphs on the same axes for comparison
     * RU: Otobrazhayet dva grafika na odnikh osyakh dlya sravneniya
     */
    void plotTwoGraphs(const QVector<double>& x1, const QVector<double>& y1,
                       const QVector<double>& x2, const QVector<double>& y2);

    /**
     * EN: Plots multiple graphs on the same axes for comparison
     * RU: Otobrazhayet neskol'ko grafikof na odnikh osyakh dlya sravneniya
     */
    void plotMultipleProfiles(const QVector<QVector<QPointF>>& profiles,
                              const QVector<QString>& labels);

    /**
     * EN: Clears all graphs from the plot
     * RU: Ochishchayet vse grafiki s polya
     */
    void clearGraphs();

    /**
     * EN: Sets axis labels
     * RU: Ustanavlivayet podpisi osey
     */
    void setAxisLabels(const QString& xLabel, const QString& yLabel);

    /**
     * EN: Sets graph title
     * RU: Ustanavlivayet zagolovok grafika
     */
    void setTitle(const QString& title);

    /**
     * EN: Enables/disables legend
     * RU: Vklyuchayet/otklyuchayet legendu
     */
    void setLegendVisible(bool visible);

    /**
     * EN: Plots multiple profiles with zero crossings marked
     * RU: Otobrazhayet neskol'ko profiliey s otmechennymi perekhodami cherez nol'
     */
    void plotMultipleProfilesWithZeroCrossings(const QVector<QVector<QPointF>>& profiles,
                                               const QVector<QString>& labels,
                                               const QVector<QVector<double>>& zeroCrossings);

    /**
     * EN: Gets pointer to QCustomPlot widget for advanced operations
     * RU: Vozvrashchaet ukazatel' na vidzhet QCustomPlot dlya rasshirennykh operatsiy
     */
    QCustomPlot* getPlot();

private:
    /**
     * EN: Updates axis ranges based on data
     * RU: Obnovlyayet diapazony osey na osnove dannykh
     */
    void updateAxisRanges(const QVector<double>& x, const QVector<double>& y);
    void updateAxisRangesForTwoGraphs(const QVector<double>& x1, const QVector<double>& y1,
                                      const QVector<double>& x2, const QVector<double>& y2);

    Ui::GraphWidget *ui;
};

#endif // GRAPHWIDGET_H