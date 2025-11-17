#ifndef IMAGESHOWCASEWIDGET_H
#define IMAGESHOWCASEWIDGET_H

#include <QWidget>
#include <QPixmap>
#include <QPainter>
#include <QMouseEvent>
#include <QPoint>
#include <QDebug>

namespace Ui {
class ImageShowcaseWidget;
}

class ImageShowcaseWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ImageShowcaseWidget(QWidget* parent = nullptr);
    explicit ImageShowcaseWidget(const QPixmap& pixmap, QWidget* parent = nullptr);
    explicit ImageShowcaseWidget(const QImage& image, QWidget* parent = nullptr);
    explicit ImageShowcaseWidget(const QString& imagePath, QWidget* parent = nullptr);

    void setImage(const QPixmap& pixmap);
    void setImage(const QImage& image);
    void setImage(const QString& imagePath);
    QPixmap getImage() const { return m_pixmap; }

    void setKeepAspectRatio(bool keep);
    void setBackgroundColor(const QColor& color);
    void setShowBorder(bool show);

    QPoint widgetToImageCoordinates(const QPoint& widgetPos) const;
    QPoint imageToWidgetCoordinates(const QPoint& imagePos) const;

signals:
    void imageClicked(const QPoint& imagePosition);
    void imageClicked(int imageX, int imageY);
    void widgetClicked(const QPoint& widgetPosition);
    void imageChanged();

protected:
    void paintEvent(QPaintEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;

private:
    void initialize();
    void updateImageRect();

    QPixmap m_pixmap;
    QRect m_imageRect;
    bool m_keepAspectRatio = true;
    QColor m_backgroundColor = Qt::lightGray;
    bool m_showBorder = true;
    Ui::ImageShowcaseWidget *ui;
};

#endif // IMAGESHOWCASEWIDGET_H
