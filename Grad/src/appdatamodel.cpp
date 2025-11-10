#include "appdatamodel.h"

AppDataModel::AppDataModel(QObject *parent)
{

}

QImage AppDataModel::currentImage() const
{
    return m_currentImage;
}

void AppDataModel::setCurrentImage(const QImage &newCurrentImage)
{
    m_originalImage = newCurrentImage;
    m_currentImage = newCurrentImage;
    emit currentImageChanged(newCurrentImage);
}

QImage AppDataModel::originalImage() const
{
    return m_originalImage;
}

double **AppDataModel::imageMatrix() const
{
    return m_imageMatrix;
}

void AppDataModel::setImageMatrix(double **newImageMatrix)
{
    m_imageMatrix = newImageMatrix;
}

double **AppDataModel::matForProfile() const
{
    return m_matForProfile;
}

void AppDataModel::setMatForProfile(double **newMatForProfile)
{
    m_matForProfile = newMatForProfile;
}

double AppDataModel::sigma() const
{
    return m_sigma;
}

void AppDataModel::setSigma(double newSigma)
{
    m_sigma = newSigma;
}

int AppDataModel::radius() const
{
    return m_radius;
}

void AppDataModel::setRadius(int newRadius)
{
    m_radius = newRadius;
}

double AppDataModel::redMax() const
{
    return m_redMax;
}

void AppDataModel::setRedMax(double newRedMax)
{
    m_redMax = newRedMax;
}

double AppDataModel::blueMax() const
{
    return m_blueMax;
}

void AppDataModel::setBlueMax(double newBlueMax)
{
    m_blueMax = newBlueMax;
}
