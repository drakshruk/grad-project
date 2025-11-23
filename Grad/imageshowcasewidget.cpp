#include "imageshowcasewidget.h"

ImageShowcaseWidget::ImageShowcaseWidget(QWidget* parent)
    : QWidget(parent)
{
    initialize();
}

ImageShowcaseWidget::ImageShowcaseWidget(const QPixmap& pixmap, QWidget* parent)
    : QWidget(parent)
    , m_pixmap(pixmap)
{
    initialize();
}

ImageShowcaseWidget::ImageShowcaseWidget(const QImage& image, QWidget* parent)
    : QWidget(parent)
    , m_pixmap(QPixmap::fromImage(image))
{
    initialize();
}

ImageShowcaseWidget::ImageShowcaseWidget(const QString& imagePath, QWidget* parent)
    : QWidget(parent)
{
    initialize();
    setImage(imagePath);
}

void ImageShowcaseWidget::initialize()
{
    setMinimumSize(100, 100);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    setFocusPolicy(Qt::StrongFocus);

    updateImageRect();
}

void ImageShowcaseWidget::setImage(const QPixmap& pixmap)
{
    m_pixmap = pixmap;
    updateImageRect();
    update();
    emit imageChanged();
}

void ImageShowcaseWidget::setImage(const QImage& image)
{
    setImage(QPixmap::fromImage(image));
}

void ImageShowcaseWidget::setImage(const QString& imagePath)
{
    QPixmap pixmap(imagePath);
    if (pixmap.isNull()) {
        qWarning() << "Failed to load image from:" << imagePath;
        m_pixmap = QPixmap();
    } else {
        m_pixmap = pixmap;
    }
    updateImageRect();
    update();
    emit imageChanged();
}

void ImageShowcaseWidget::setKeepAspectRatio(bool keep)
{
    if (m_keepAspectRatio != keep) {
        m_keepAspectRatio = keep;
        updateImageRect();
        update();
    }
}

void ImageShowcaseWidget::setBackgroundColor(const QColor& color)
{
    if (m_backgroundColor != color) {
        m_backgroundColor = color;
        update();
    }
}

void ImageShowcaseWidget::setShowBorder(bool show)
{
    if (m_showBorder != show) {
        m_showBorder = show;
        update();
    }
}

QPoint ImageShowcaseWidget::widgetToImageCoordinates(const QPoint& widgetPos) const
{
    if (m_pixmap.isNull() || m_imageRect.isEmpty())
        return QPoint(-1, -1);

    if (m_imageRect.contains(widgetPos)) {
        double scaleX = static_cast<double>(m_pixmap.width()) / m_imageRect.width();
        double scaleY = static_cast<double>(m_pixmap.height()) / m_imageRect.height();

        int imgX = static_cast<int>((widgetPos.x() - m_imageRect.x()) * scaleX);
        int imgY = static_cast<int>((widgetPos.y() - m_imageRect.y()) * scaleY);

        imgX = qBound(0, imgX, m_pixmap.width() - 1);
        imgY = qBound(0, imgY, m_pixmap.height() - 1);

        return QPoint(imgX, imgY);
    }

    return QPoint(-1, -1);
}

QPoint ImageShowcaseWidget::imageToWidgetCoordinates(const QPoint& imagePos) const
{
    if (m_pixmap.isNull() || m_imageRect.isEmpty())
        return QPoint(-1, -1);

    if (imagePos.x() >= 0 && imagePos.x() < m_pixmap.width() &&
        imagePos.y() >= 0 && imagePos.y() < m_pixmap.height()) {

        double scaleX = static_cast<double>(m_imageRect.width()) / m_pixmap.width();
        double scaleY = static_cast<double>(m_imageRect.height()) / m_pixmap.height();

        int widgetX = m_imageRect.x() + static_cast<int>(imagePos.x() * scaleX);
        int widgetY = m_imageRect.y() + static_cast<int>(imagePos.y() * scaleY);

        return QPoint(widgetX, widgetY);
    }

    return QPoint(-1, -1);
}

void ImageShowcaseWidget::updateImageRect()
{
    if (m_pixmap.isNull()) {
        m_imageRect = QRect();
        return;
    }

    QRect widgetRect = rect();

    if (m_keepAspectRatio) {
        QSize pixmapSize = m_pixmap.size();
        pixmapSize.scale(widgetRect.size(), Qt::KeepAspectRatio);

        m_imageRect = QRect(0, 0, pixmapSize.width(), pixmapSize.height());
        m_imageRect.moveCenter(widgetRect.center());
    } else {
        m_imageRect = widgetRect;
    }
}

void ImageShowcaseWidget::paintEvent(QPaintEvent* event)
{
    Q_UNUSED(event)

    QPainter painter(this);

    painter.fillRect(rect(), m_backgroundColor);

    if (m_showBorder) {
        painter.setPen(QPen(Qt::darkGray, 2));
        painter.drawRect(rect().adjusted(1, 1, -1, -1));
    }

    if (!m_pixmap.isNull()) {
        painter.drawPixmap(m_imageRect, m_pixmap);

        painter.setPen(QPen(Qt::black, 1, Qt::DotLine));
        painter.drawRect(m_imageRect);
    } else {
        painter.setPen(Qt::darkGray);
        painter.drawText(rect(), Qt::AlignCenter, "No Image Loaded");
    }
}

void ImageShowcaseWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        QPoint widgetPos = event->pos();
        QPoint imagePos = widgetToImageCoordinates(widgetPos);

        qDebug() << "Widget clicked at:" << widgetPos;
        qDebug() << "Image coordinates:" << imagePos;

        emit widgetClicked(widgetPos);

        if (imagePos.x() >= 0 && imagePos.y() >= 0) {
            emit imageClicked(imagePos);
            emit imageClicked(imagePos.x(), imagePos.y());
        }
    }

    QWidget::mousePressEvent(event);
}

void ImageShowcaseWidget::resizeEvent(QResizeEvent* event)
{
    QWidget::resizeEvent(event);
    updateImageRect();
    update();
}
