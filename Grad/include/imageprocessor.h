#ifndef IMAGEPROCESSOR_H
#define IMAGEPROCESSOR_H

#include <imageprocessingtypes.h>

namespace ImageProcessor {

    const QImage substractTwoImages(const QImage &Im1, const QImage &Im2);

    const QImage addTwoImages(const QImage &Im1, const QImage &Im2);

    const QImage multiplyTwoImages(const QImage &Im1, const QImage &Im2);

    const QImage divideTwoImages(const QImage &Im1, const QImage &Im2);

    const QImage minTwoImages(const QImage &Im1, const QImage &Im2);

    const QImage maxTwoImages(const QImage &Im1, const QImage &Im2);

    const QImage AND_TwoImages(const QImage &Im1, const QImage &Im2);

    const QImage OR_TwoImages(const QImage &Im1, const QImage &Im2);

    const QImage XOR_TwoImagess(const QImage &Im1, const QImage &Im2);

    const QImage averageTwoImages(const QImage &Im1, const QImage &Im2);

    const QImage differenceTwoImages(const QImage &Im1, const QImage &Im2);

    const QImage transparentZeroTwoImages(const QImage &Im1, const QImage &Im2);

    bool containsVP(QVector<QPoint> &vec, QPoint p);

    void selectEdge(QPoint pos, Matrix2D<int> &iiImAttMat, QVector<QPoint> &selectedEdge);

//    const Matrix2D<double> substractTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> addTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> multiplyTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> divideTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> minTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> maxTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> AND_TwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> OR_TwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> XOR_TwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> averageTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> differenceTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

//    const Matrix2D<double> transparentZeroTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2);

    template <typename T, typename Operation>
    const Matrix2D<T> elementWiseOperation(const Matrix2D<T> &mat1, const Matrix2D<T> &mat2, Operation op)
    {
        size_t w = min(mat1.size(), mat2.size());
        size_t h = min(mat1[0].size(), mat2[0].size());
        Matrix2D<T> res(w, std::vector<T>(h, 0));
        for(size_t i = 0; i < w; i++) {
            const auto& row1 = mat1[i];
            const auto& row2 = mat2[i];
            auto& resRow = res[i];

            for(size_t j = 0; j < h; j++) {
                resRow[j] = op(row1[j], row2[j]);
            }
        }
        return res;
    }

    const ProfileResult buildProfile001(const ProfileParameters &params, const Matrix2D<double> &ddImMat);

    const ProfileResult buildProfile002(const ProfileParameters &params, const Matrix2D<double> &ddImMat, int convRad);

    const ProfileResult cutProfile(const ProfileResult &profile, int step);

    const ProfileResult cutProfile(const ProfileResult &profile, int left, int right);

    void alignProfiles(QVector<QPointF> &profile1, QVector<QPointF> &profile2, int len);

    double differenceOfTwoProfiles(const ProfileResult &profile1, const ProfileResult &profile2);

    const ProfileResult reInterpolateProfile(const ProfileResult &profile, int newLen);

    const Matrix2D<double> fromGrayImage(const QImage &image);

    const QImage toGrayImage(const Matrix2D<double> &ddImageMat);

    const QImage toBlueRedImage(const Matrix2D<double> &ddImageMat, double dRedMax, double dBlueMax);

    const QImage combineImageWithEdge(const Matrix2D<double> &ddImageMat, const Matrix2D<int> &iiImageAttMat);

    const QImage combineImageWithEdge(const QImage &Image, const Matrix2D<int> &iiImageAttMat);

    const QImage combineImageWithProfile(const QImage &im1, const QVector<QPoint> &profileIndices);

    const QImage convImage(const QImage &image, const Matrix2D<double> &ddConvCore);

    const Matrix2D<double> convMat(const Matrix2D<double> &ddMat, const Matrix2D<double> &ddConvMat);

    const Matrix2D<double> getGauss(int iXsize, int iYsize, double dSigma);

    const Matrix2D<double> getXGradCore(int iXsize, int iYsize, double dSigma);

    const Matrix2D<double> getYGradCore(int iXsize, int iYsize, double dSigma);

    const Matrix2D<double> getLapl(int iXsize, int iYsize, double dSigma);

    const Matrix2D<double> findEdges(const Matrix2D<double> &ddImageMat, int iRad, Matrix2D<int> &iiImageAttMat);

    const Matrix2D<double> findEdges(const Matrix2D<double> &ddImageMat, int iRad);

    const Matrix2D<double> findEdgesRB(const Matrix2D<double> &ddImageMat, int iRad, Matrix2D<int> &iiImageAttMat);

    const QImage gaussianEdgeDetection(const Matrix2D<double> &ddImageMat, double dSigma, int iRad, Matrix2D<int> &iiImageAttMat, Matrix2D<double> &ddMatForProfile);

    const QImage gaussianEdgeDetection(const Matrix2D<double> &ddImageMat, double dSigma, int iRad);

    const QImage sampleTwoHollows();

    const QImage sampleTwoHollowsBig();
};

class ImageData {
public:
    QImage currentImage;
    QImage originalImage; // If you need both
    Matrix2D<double> imageMatrix;
    Matrix2D<int> attributeMatrix;
    Matrix2D<double> testMatrix;

    ImageData(const QImage& image){
        currentImage = image;
        originalImage = image;

        imageMatrix = ImageProcessor::fromGrayImage(image);

        const int w = image.width();
        const int h = image.height();

        attributeMatrix = Matrix2D<int>(w, std::vector<int>(h, 0));
        testMatrix = Matrix2D<double>(w, std::vector<double>(h, 0.0));
    }

    void updateOriginalFromImage(const QImage& newImage) {
        currentImage = newImage;
        originalImage = newImage;

        imageMatrix = ImageProcessor::fromGrayImage(newImage);

        // Initialize other matrices
        const int w = newImage.width();
        const int h = newImage.height();

        attributeMatrix = Matrix2D<int>(w, std::vector<int>(h, 0));
        testMatrix = Matrix2D<double>(w, std::vector<double>(h, 0.0));
    }

    void updateCurrentFromImage(const QImage& newImage) {
        currentImage = newImage;

        imageMatrix = ImageProcessor::fromGrayImage(newImage);

        // Initialize other matrices
        const int w = newImage.width();
        const int h = newImage.height();

        attributeMatrix = Matrix2D<int>(w, std::vector<int>(h, 0));
        testMatrix = Matrix2D<double>(w, std::vector<double>(h, 0.0));
    }

    void clear() {
        imageMatrix.clear();
        attributeMatrix.clear();
        testMatrix.clear();
        currentImage = QImage();
        originalImage = QImage();
    }
};

#endif // IMAGEPROCESSOR_H
