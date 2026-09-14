#ifndef IMAGEPROCESSOR_H
#define IMAGEPROCESSOR_H

#include "imageprocessingtypes.h"

namespace ImageProcessor {

    // ============================================================================
    // Utility Functions / Vspomogatel'nye funkcii
    // ============================================================================

    /**
     * @brief Checks if vector contains point p
     * @param vec Vector of points to search
     * @param p Point to find
     * @return true if point exists, false otherwise
     */
    bool containsVP(const QVector<QPoint>& vec, const QPoint& p);

    /**
     * @brief Selects connected edge pixels starting from position pos
     * @param pos Starting position (must be an edge pixel)
     * @param iiImAttMat Attribute matrix (modified)
     * @param selectedEdge Output vector of selected edge points
     */
    void selectEdge(const QPoint& pos, Matrix2D<int>& iiImAttMat, QVector<QPoint>& selectedEdge);

    // ============================================================================
    // Profile Building Functions / Funkcii postroeniya profilya
    // ============================================================================

    /**
     * @brief Builds profile perpendicular to edge by sampling along gradient direction
     * @param params Profile parameters (position, gradient, sigma)
     * @param ddImMat Input image matrix
     * @return ProfileResult containing sampled points and indices
     */
    ProfileResult buildProfile001(const ProfileParameters& params, const Matrix2D<double>& ddImMat);

    /**
     * @brief Builds profile with Laplacian convolution
     * @param params Profile parameters
     * @param ddImMat Input image matrix
     * @param convRad Convolution kernel radius
     * @return ProfileResult with convolved values
     */
    ProfileResult buildProfile002(const ProfileParameters& params, const Matrix2D<double>& ddImMat, int convRad);

    /**
     * @brief Cuts profile by stepping 'step' pixels from center to local min/max
     * @param profile Input profile
     * @param step Number of pixels to step beyond local extrema
     * @return Cut profile
     * @note Currently unused - kept for future use
     */
    QVector<QPointF> cutProfile(const QVector<QPointF>& profile, int step);

    /**
     * @brief Cuts profile between left and right indices
     * @param profile Input profile
     * @param left Left boundary index
     * @param right Right boundary index
     * @return Cut profile
     */
    QVector<QPointF> cutProfile(const QVector<QPointF>& profile, int left, int right);

    /**
     * @brief Calculates residual between two profiles (RMS error)
     * @param profile1 First profile
     * @param profile2 Second profile
     * @return Root mean square difference
     */
    double calculateResidualOfProfiles(const QVector<QPointF>& profile1, const QVector<QPointF>& profile2);

    /**
     * @brief Calculates sum of squared differences between profiles
     * @param profile1 First profile
     * @param profile2 Second profile
     * @return Sum of squared differences
     * @note Currently unused - kept for future use
     */
    double differenceOfTwoProfiles(const QVector<QPointF>& profile1, const QVector<QPointF>& profile2);

    /**
     * @brief Reinterpolates profile to new length using linear interpolation
     * @param profile Input profile
     * @param newLen Desired output length
     * @return Reinterpolated profile
     */
    QVector<QPointF> reInterpolateProfile(const QVector<QPointF>& profile, int newLen);

    /**
     * @brief Reinterpolates profile from leftBound to rightBound to new length
     * @param profile Input profile
     * @param leftBound Left boundary (0 to size-1)
     * @param rightBound Right boundary (0 to size-1)
     * @param newLen Desired output length
     * @return Reinterpolated profile
     */
    QVector<QPointF> reInterpolateProfile(const QVector<QPointF>& profile, double leftBound, double rightBound, int newLen);

    // ============================================================================
    // Image Conversion Functions / Funkcii preobrazovaniya izobrazheniy
    // ============================================================================

    /**
     * @brief Converts QImage to grayscale double matrix
     * @param image Input image
     * @return 2D matrix of double values (0-255)
     */
    Matrix2D<double> fromGrayImage(const QImage& image);

    /**
     * @brief Converts double matrix to grayscale QImage
     * @param ddImageMat Input matrix
     * @return Grayscale QImage
     */
    QImage toGrayImage(const Matrix2D<double>& ddImageMat);

    /**
     * @brief Converts matrix to red/blue image (positive=red, negative=blue)
     * @param ddImageMat Input matrix
     * @param dRedMax Maximum red intensity scaling
     * @param dBlueMax Maximum blue intensity scaling
     * @return RGB image with red/blue pixels
     */
    QImage toBlueRedImage(const Matrix2D<double>& ddImageMat, double dRedMax, double dBlueMax);

    /**
     * @brief Combines image with edge overlay
     * @param ddImageMat Image matrix
     * @param iiImageAttMat Attribute matrix (edges marked)
     * @return Combined image
     */
    QImage combineImageWithEdge(const Matrix2D<double>& ddImageMat, const Matrix2D<int>& iiImageAttMat);

    /**
     * @brief Combines QImage with edge overlay
     * @param Image Input image
     * @param iiImageAttMat Attribute matrix
     * @return Combined image
     */
    QImage combineImageWithEdge(const QImage& Image, const Matrix2D<int>& iiImageAttMat);

    /**
     * @brief Combines image with profile indices (colors them cyan)
     * @param im1 Input image
     * @param profileIndices Points to highlight
     * @return Combined image
     * @note Currently unused - kept for future use
     */
    QImage combineImageWithProfile(const QImage& im1, const QVector<QPoint>& profileIndices);

    // ============================================================================
    // Convolution and Kernel Functions / Funkcii svertki i yader
    // ============================================================================

    /**
     * @brief Convolves QImage with kernel
     * @param image Input image
     * @param ddConvCore Convolution kernel
     * @return Convolved image
     */
    QImage convImage(const QImage& image, const Matrix2D<double>& ddConvCore);

    /**
     * @brief Convolves double matrix with kernel
     * @param ddMat Input matrix
     * @param ddConvMat Convolution kernel
     * @return Convolved matrix
     */
    Matrix2D<double> convMat(const Matrix2D<double>& ddMat, const Matrix2D<double>& ddConvMat);

    /**
     * @brief Generates Gaussian kernel
     * @param iXsize Kernel width
     * @param iYsize Kernel height
     * @param dSigma Gaussian sigma
     * @return Gaussian kernel matrix
     */
    Matrix2D<double> getGauss(int iXsize, int iYsize, double dSigma);

    /**
     * @brief Generates X-gradient kernel (Gaussian derivative)
     * @param iXsize Kernel width
     * @param iYsize Kernel height
     * @param dSigma Gaussian sigma
     * @return X-gradient kernel
     */
    Matrix2D<double> getXGradCore(int iXsize, int iYsize, double dSigma);

    /**
     * @brief Generates Y-gradient kernel (Gaussian derivative)
     * @param iXsize Kernel width
     * @param iYsize Kernel height
     * @param dSigma Gaussian sigma
     * @return Y-gradient kernel
     */
    Matrix2D<double> getYGradCore(int iXsize, int iYsize, double dSigma);

    /**
     * @brief Generates Laplacian of Gaussian kernel
     * @param iXsize Kernel width
     * @param iYsize Kernel height
     * @param dSigma Gaussian sigma
     * @return LoG kernel
     */
    Matrix2D<double> getLapl(int iXsize, int iYsize, double dSigma);

    // ============================================================================
    // Edge Detection Functions / Funkcii detekcii granic
    // ============================================================================

    /**
     * @brief Finds edge pixels (where positive and negative neighbors exist)
     * @param ddImageMat Input matrix
     * @param iRad Search radius
     * @param iiImageAttMat Output attribute matrix (modified)
     * @return Edge magnitude matrix
     */
    Matrix2D<double> findEdges(const Matrix2D<double>& ddImageMat, int iRad, Matrix2D<int>& iiImageAttMat);

    /**
     * @brief Finds edge pixels without attribute matrix
     * @param ddImageMat Input matrix
     * @param iRad Search radius
     * @return Edge magnitude matrix
     */
    Matrix2D<double> findEdges(const Matrix2D<double>& ddImageMat, int iRad);

    /**
     * @brief Finds edges with red/blue encoding (positive/negative)
     * @param ddImageMat Input matrix
     * @param iRad Search radius
     * @param iiImageAttMat Output attribute matrix
     * @return Edge matrix with signed values (-1, 0, 1)
     * @note Currently unused - kept for future use
     */
    Matrix2D<double> findEdgesRB(const Matrix2D<double>& ddImageMat, int iRad, Matrix2D<int>& iiImageAttMat);

    /**
     * @brief Gaussian edge detection with attribute matrix output
     * @param ddImageMat Input image matrix
     * @param dSigma Gaussian sigma
     * @param iRad Kernel radius
     * @param iiImageAttMat Output attribute matrix (modified)
     * @return Edge detection result image
     */
    QImage gaussianEdgeDetection(const Matrix2D<double>& ddImageMat, double dSigma, int iRad,
                                  Matrix2D<int>& iiImageAttMat);

    /**
     * @brief Simplified Gaussian edge detection
     * @param ddImageMat Input image matrix
     * @param dSigma Gaussian sigma
     * @param iRad Kernel radius
     * @return Edge detection result image
     */
    QImage gaussianEdgeDetection(const Matrix2D<double>& ddImageMat, double dSigma, int iRad);

    // ============================================================================
    // Sample Image Generation / Generaciya testovyh izobrazheniy
    // ============================================================================

    /**
     * @brief Creates test image with two hollow circles (200x200)
     * @return Sample image
     */
    QImage sampleTwoHollows();

    /**
     * @brief Creates test image with two hollow circles (400x550)
     * @return Sample image
     */
    QImage sampleTwoHollowsBig();

    // ============================================================================
    // Edge Refinement / Utochnenie granic
    // ============================================================================

    /**
     * @brief Refines a single edge point using PCR algorithm
     * @param n0 Initial X coordinate
     * @param m0 Initial Y coordinate
     * @param params Refinement parameters
     * @return RefinementResult with refined position and statistics
     */
    RefinementResult refineSinglePoint001(int n0, int m0, const RefinementParameters& params);
    RefinementResult refineSinglePoint002(int n0, int m0, const RefinementParameters& params);



    // Вычисление остатка между двумя профилями (для градиентного спуска)
    double computeResidual(const QVector<double>& yy1,
                                           const QVector<QPointF>& yP2_full,
                                           int XL2, int XR2, int NN);

    double interpolateBilinear(const Matrix2D<double>& img, double x, double y);

    // ============================================================================
    // Element-wise Operations / Poelementnye operacii
    // ============================================================================

    /**
     * @brief Applies element-wise operation to two matrices
     * @param mat1 First matrix
     * @param mat2 Second matrix
     * @param op Operation functor
     * @return Result matrix
     */
    template <typename T, typename Operation>
    Matrix2D<T> elementWiseOperation(const Matrix2D<T>& mat1, const Matrix2D<T>& mat2, Operation op)
    {
        if (mat1.empty() || mat2.empty() || mat1[0].empty() || mat2[0].empty()) {
            return Matrix2D<T>();
        }

        size_t rows = std::min(mat1.size(), mat2.size());
        size_t cols = std::min(mat1[0].size(), mat2[0].size());

        Matrix2D<T> res(rows, std::vector<T>(cols, T(0)));

        for(size_t i = 0; i < rows; i++) {
            const auto& row1 = mat1[i];
            const auto& row2 = mat2[i];
            auto& resRow = res[i];

            for(size_t j = 0; j < cols; j++) {
                resRow[j] = op(row1[j], row2[j]);
            }
        }
        return res;
    }
    // ============================================================================
    // Profile Building Functions (New) / Funkcii postroeniya profilya (Novye)
    // ============================================================================

    /**
     * EN: Builds profile along a line between two points using Bresenham's algorithm
     * RU: Stroit profil' vdol' linii mezhdu dvumya tochkami s ispol'zovaniyem algoritma Brezenhema
     * @param p1 First point / Pervaya tochka
     * @param p2 Second point / Vtoraya tochka
     * @param image Input image / Vkhodnoye izobrazheniye
     * @return Vector of intensity values along the line / Vektor znacheniy intensivnosti vdol' linii
     */
    QVector<double> buildProfileBetweenPoints(const QPoint& p1, const QPoint& p2, const QImage& image);

    /**
     * EN: Builds profile along a line between two points with interpolation for subpixel accuracy
     * RU: Stroit profil' vdol' linii mezhdu dvumya tochkami s interpolyatsiyey dlya subpiksel'noy tochnosti
     * @param p1 First point (floating point) / Pervaya tochka (s plavayushchey tochkoy)
     * @param p2 Second point (floating point) / Vtoraya tochka (s plavayushchey tochkoy)
     * @param image Input image / Vkhodnoye izobrazheniye
     * @param numSamples Number of samples along the line / Kolichestvo otschetov vdol' linii
     * @return Vector of intensity values along the line / Vektor znacheniy intensivnosti vdol' linii
     */
    QVector<double> buildProfileBetweenPoints(const QPointF& p1, const QPointF& p2,
                                              const QImage& image, int numSamples = 100);

    // ============================================================================
    // Image Statistics Functions (New) / Funkcii statistiki izobrazheniy (Novye)
    // ============================================================================

    /**
     * EN: Structure containing comprehensive image statistics
     * RU: Struktura, soderzhashchaya vsestoronyuyu statistiku izobrazheniya
     */
    struct ImageStatistics {
        double mean = 0.0;           // Mean intensity / Srednyaya intensivnost'
        double median = 0.0;         // Median intensity / Mediana intensivnosti
        double min = 0.0;            // Minimum intensity / Minimal'naya intensivnost'
        double max = 0.0;            // Maximum intensity / Maksimal'naya intensivnost'
        double sum = 0.0;            // Sum of all intensities / Summa vsekh intensivnostey
        double variance = 0.0;       // Variance / Dispersiya
        double stdDev = 0.0;         // Standard deviation / Srednekvadraticheskoye otkloneniye
        double skewness = 0.0;       // Skewness / Asimmetriya
        double kurtosis = 0.0;       // Kurtosis / Ekstsess
        double entropy = 0.0;        // Entropy / Entropiya
        long long pixelCount = 0;    // Total number of pixels / Obshcheye kolichestvo pikseley

        // Histogram data (256 bins for grayscale) / Dannyye gistogrammy (256 binov dlya seroy shkaly)
        QVector<long long> histogram;  // 256 bins / 256 binov

        // Additional stats for edge detection / Dopolnitel'naya statistika dlya detektsii granits
        double edgeRatio = 0.0;       // Ratio of edge pixels to total pixels / Otnosheniye granichnykh pikseley k obshchemu kolichestvu
        int edgePixelCount = 0;       // Number of edge pixels / Kolichestvo granichnykh pikseley

        // Format as string for display / Formatirovat' kak stroku dlya otobrazheniya
        QString toString() const;
    };

    /**
     * EN: Computes comprehensive statistics for an image
     * RU: Vychislyayet vsestoronyuyu statistiku dlya izobrazheniya
     * @param image Input image / Vkhodnoye izobrazheniye
     * @param attributeMatrix Optional attribute matrix for edge statistics (can be empty)
     * @return ImageStatistics structure with all computed values
     */
    ImageStatistics computeImageStatistics(const QImage& image,
                                           const Matrix2D<int>& attributeMatrix = Matrix2D<int>());

    /**
     * EN: Computes statistics for a double matrix (for gradient images, etc.)
     * RU: Vychislyayet statistiku dlya matritsy double (dlya gradientnykh izobrazheniy, itd.)
     * @param matrix Input matrix / Vkhodnaya matritsa
     * @return ImageStatistics structure
     */
    ImageStatistics computeMatrixStatistics(const Matrix2D<double>& matrix);

    /**
     * EN: Saves statistics to a text file
     * RU: Sohranyayet statistiku v tekstovyy fayl
     * @param stats Statistics to save / Statistika dlya sohraneniya
     * @param filename Output file path / Put' k vykhodnomu faylu
     * @return true if successful, false otherwise
     */
    bool saveStatisticsToFile(const ImageStatistics& stats, const QString& filename);

}
// namespace ImageProcessor

// ============================================================================
// ImageData Class - Data container for image and associated matrices
// Note: Consider using this class instead of storing matrices separately in MainWindow
// ============================================================================

class ImageData {
public:
    QImage currentImage;
    QImage originalImage;
    Matrix2D<double> imageMatrix;
    Matrix2D<int> attributeMatrix;
    Matrix2D<double> testMatrix;  // For profile testing

    // Constructor
    explicit ImageData(const QImage& image) {
        updateFromImage(image);
    }

    // Default constructor
    ImageData() = default;

    // Copy constructor
    ImageData(const ImageData& other) = default;

    // Move constructor
    ImageData(ImageData&& other) noexcept = default;

    // Copy assignment
    ImageData& operator=(const ImageData& other) = default;

    // Move assignment
    ImageData& operator=(ImageData&& other) noexcept = default;

    // Destructor - matrices will be cleaned up automatically (std::vector handles it)
    ~ImageData() = default;

    /**
     * @brief Updates all data from a new image
     * @param newImage New image to load
     */
    void updateFromImage(const QImage& newImage) {
        currentImage = newImage;
        originalImage = newImage;

        imageMatrix = ImageProcessor::fromGrayImage(newImage);

        const int w = newImage.width();
        const int h = newImage.height();

        // Reinitialize matrices with correct dimensions
        attributeMatrix = Matrix2D<int>(w, std::vector<int>(h, 0));
        testMatrix = Matrix2D<double>(w, std::vector<double>(h, 0.0));
    }

    /**
     * @brief Updates only current image (keeps original)
     * @param newImage New current image
     */
    void updateCurrentFromImage(const QImage& newImage) {
        currentImage = newImage;
        imageMatrix = ImageProcessor::fromGrayImage(newImage);

        const int w = newImage.width();
        const int h = newImage.height();

        attributeMatrix = Matrix2D<int>(w, std::vector<int>(h, 0));
        testMatrix = Matrix2D<double>(w, std::vector<double>(h, 0.0));
    }

    /**
     * @brief Clears all data
     */
    void clear() {
        imageMatrix.clear();
        attributeMatrix.clear();
        testMatrix.clear();
        currentImage = QImage();
        originalImage = QImage();
    }

    /**
     * @brief Checks if data is valid (has an image loaded)
     * @return true if image is valid
     */
    bool isValid() const {
        return !currentImage.isNull() && !imageMatrix.empty();
    }

    /**
     * @brief Gets image dimensions
     */
    QSize size() const {
        return currentImage.size();
    }

    int width() const { return currentImage.width(); }
    int height() const { return currentImage.height(); }
};

#endif // IMAGEPROCESSOR_H
