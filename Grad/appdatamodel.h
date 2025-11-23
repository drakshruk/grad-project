#ifndef APPDATAMODEL_H
#define APPDATAMODEL_H

#include <QImage>
#include <QObject>

class AppDataModel : public QObject
{
    Q_OBJECT

public:

    explicit AppDataModel(QObject *parent = nullptr);

    QImage currentImage() const;
    void setCurrentImage(const QImage &newCurrentImage);


    QImage originalImage() const;
    void setOriginalImage(const QImage &newOriginalImage);

    QVector<QPoint> selectedEdge() const;
    void setSelectedEdge(const QVector<QPoint> &newSelectedEdge);

    double** imageMatrix() const;
    void setImageMatrix(double** newImageMatrix);

    double** matForProfile() const;
    void setMatForProfile(double** newMatForProfile);

    double sigma() const;
    void setSigma(double newSigma);

    int radius() const;
    void setRadius(int newRadius);

    double redMax() const;
    void setRedMax(double newRedMax);

    double blueMax() const;
    void setBlueMax(double newBlueMax);


signals:
    void currentImageChanged(const QImage &image);
//    void sigmaChanged(double sigma);
//    void radiusChanged(int radius);
//    void redMaxChanged(double redMax);
//    void blueMaxChanged(double blueMax);

private:
    QImage m_currentImage;
    QImage m_originalImage;
    QVector<QPoint> m_selectedEdge;
    double** m_imageMatrix = nullptr;
    int** m_attribiteMatrix = nullptr;
    double** m_matForProfile = nullptr;

    // Constants/Parameters
    double m_sigma = 5.0;
    int m_radius = 10;
    double m_redMax = 255.0;
    double m_blueMax = 255.0;
};

#endif // APPDATAMODEL_H
