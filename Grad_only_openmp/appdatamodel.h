#ifndef APPDATAMODEL_H
#define APPDATAMODEL_H

#include <QImage>
#include <QObject>
#include <QVector>
#include <QPoint>
#include "imageprocessingtypes.h"

/**
 * @brief Application Data Model - Central data storage for the application
 * @brief Model' dannykh prilozheniya - Tsentral'noye khranilishche dannykh
 *
 * EN: This class serves as the single source of truth for all application data.
 *     It uses Qt's signal/slot mechanism to notify views when data changes.
 * RU: Etot klass sluzhit yedinstvennym istochnikom istiny dlya vsekh dannykh prilozheniya.
 *     Ispol'zuyet mekhanizm signalov/slotov Qt dlya opoveshcheniya vidov ob izmenenii dannykh.
 */
class AppDataModel : public QObject
{
    Q_OBJECT

public:
    // Constructor / Destructor / Konstruktor / Destruktor
    explicit AppDataModel(QObject *parent = nullptr);
    ~AppDataModel();

    // Disable copy to prevent accidental duplication (enable if needed)
    // Zapreshchayem kopirovaniye dlya predotvrashcheniya sluchaynogo duplirovaniya
    AppDataModel(const AppDataModel&) = delete;
    AppDataModel& operator=(const AppDataModel&) = delete;

    // Enable move semantics / Razreshayem peremeshcheniye
    AppDataModel(AppDataModel&&) noexcept = default;
    AppDataModel& operator=(AppDataModel&&) noexcept = default;

    // ========================================================================
    // Image Data / Dannyye izobrazheniya
    // ========================================================================

    /**
     * @brief Gets current working image / Vozvrashchayet tekushcheye rabocheye izobrazheniye
     */
    QImage currentImage() const { return m_currentImage; }

    /**
     * @brief Sets current working image and notifies observers
     * @brief Ustanavlivayet tekushcheye rabocheye izobrazheniye i opoveshchayet nablyudateley
     */
    void setCurrentImage(const QImage &newCurrentImage);

    /**
     * @brief Gets original (unmodified) image / Vozvrashchayet original'noye (neizmenennoye) izobrazheniye
     */
    QImage originalImage() const { return m_originalImage; }

    /**
     * @brief Sets original image / Ustanavlivayet original'noye izobrazheniye
     */
    void setOriginalImage(const QImage &newOriginalImage);

    /**
     * @brief Checks if model has valid image data / Proveryayet, soderzhit li model' korrektnyye dannyye izobrazheniya
     */
    bool hasImage() const { return !m_currentImage.isNull(); }

    /**
     * @brief Gets image dimensions / Vozvrashchayet razmery izobrazheniya
     */
    QSize imageSize() const { return m_currentImage.size(); }
    int width() const { return m_currentImage.width(); }
    int height() const { return m_currentImage.height(); }

    // ========================================================================
    // Matrix Data / Matrichnyye dannyye
    // ========================================================================

    /**
     * @brief Gets grayscale image matrix / Vozvrashchayet matritsu izobrazheniya v ottenkakh serogo
     */
    Matrix2D<double> imageMatrix() const { return m_imageMatrix; }

    /**
     * @brief Sets grayscale image matrix / Ustanavlivayet matritsu izobrazheniya v ottenkakh serogo
     */
    void setImageMatrix(const Matrix2D<double>& newImageMatrix);

    /**
     * @brief Gets attribute matrix (edge flags) / Vozvrashchayet atributnuyu matritsu (flagi granits)
     */
    Matrix2D<int> attributeMatrix() const { return m_attributeMatrix; }

    /**
     * @brief Sets attribute matrix / Ustanavlivayet atributnuyu matritsu
     */
    void setAttributeMatrix(const Matrix2D<int>& newAttributeMatrix);

    /**
     * @brief Gets test matrix (for profile calculations) / Vozvrashchayet testovuyu matritsu (dlya raschetov profilya)
     */
    Matrix2D<double> testMatrix() const { return m_testMatrix; }

    /**
     * @brief Sets test matrix / Ustanavlivayet testovuyu matritsu
     */
    void setTestMatrix(const Matrix2D<double>& newTestMatrix);

    // ========================================================================
    // Edge Selection Data / Dannyye vybora granits
    // ========================================================================

    /**
     * @brief Gets currently selected edge points / Vozvrashchayet vybrannyye tochki granits
     */
    QVector<QPoint> selectedEdge() const { return m_selectedEdge; }

    /**
     * @brief Sets selected edge points / Ustanavlivayet vybrannyye tochki granits
     */
    void setSelectedEdge(const QVector<QPoint>& newSelectedEdge);

    /**
     * @brief Clears edge selection / Ochishchayet vybor granits
     */
    void clearEdgeSelection();

    /**
     * @brief Checks if edge selection mode is active / Proveryayet, aktivirovan li rezhim vybora granits
     */
    bool isEdgeSelectionMode() const { return m_edgeSelectionMode; }

    /**
     * @brief Sets edge selection mode / Ustanavlivayet rezhim vybora granits
     */
    void setEdgeSelectionMode(bool enabled);

    // ========================================================================
    // Algorithm Parameters / Parametry algoritmov
    // ========================================================================

    /**
     * @brief Gets Gaussian sigma value / Vozvrashchayet znacheniye sigma Gaussa
     */
    double sigma() const { return m_sigma; }

    /**
     * @brief Sets Gaussian sigma value and notifies observers
     * @brief Ustanavlivayet znacheniye sigma Gaussa i opoveshchayet nablyudateley
     */
    void setSigma(double newSigma);

    /**
     * @brief Gets kernel radius / Vozvrashchayet radius yadra
     */
    int radius() const { return m_radius; }

    /**
     * @brief Sets kernel radius and notifies observers
     * @brief Ustanavlivayet radius yadra i opoveshchayet nablyudateley
     */
    void setRadius(int newRadius);

    /**
     * @brief Gets maximum red channel value for visualization / Vozvrashchayet maksimal'noye znacheniye krasnogo kanala dlya vizualizatsii
     */
    double redMax() const { return m_redMax; }

    /**
     * @brief Sets maximum red channel value / Ustanavlivayet maksimal'noye znacheniye krasnogo kanala
     */
    void setRedMax(double newRedMax);

    /**
     * @brief Gets maximum blue channel value for visualization / Vozvrashchayet maksimal'noye znacheniye sinego kanala dlya vizualizatsii
     */
    double blueMax() const { return m_blueMax; }

    /**
     * @brief Sets maximum blue channel value / Ustanavlivayet maksimal'noye znacheniye sinego kanala
     */
    void setBlueMax(double newBlueMax);

    // ========================================================================
    // Utility Methods / Vspomogatel'nyye metody
    // ========================================================================

    /**
     * @brief Updates all matrices from current image / Obnovlyayet vse matritsy iz tekushchego izobrazheniya
     */
    void updateMatricesFromCurrentImage();

    /**
     * @brief Clears all data / Ochishchayet vse dannyye
     */
    void clear();

    /**
     * @brief Resets to original image / Vosstanavlivayet original'noye izobrazheniye
     */
    void resetToOriginal();


    /**
     * EN: Validates that sigma is not larger than image_size/8
     * RU: Proveryayet, chto sigma ne bol'she chem razmer_izobrazheniya/8
     * @param sigma Value to validate / Znacheniye dlya proverki
     * @return Clamped valid sigma value / Ogranichennoye korrektnoye znacheniye sigma
     */
    double validateSigma(double sigma) const;

    /**
     * EN: Gets maximum allowed sigma based on current image size
     * RU: Poluchayet maksimal'no dopustimoye znacheniye sigma na osnove tekushchego razmera izobrazheniya
     * @return Maximum allowed sigma / Maksimal'naya dopustimaya sigma
     */
    double maxAllowedSigma() const;

signals:
    // ========================================================================
    // Data Change Signals / Signaly ob izmenenii dannykh
    // ========================================================================

    /**
     * @brief Emitted when current image changes / Ispuskayetsya pri izmenenii tekushchego izobrazheniya
     */
    void currentImageChanged(const QImage &image);

    /**
     * @brief Emitted when original image changes / Ispuskayetsya pri izmenenii original'nogo izobrazheniya
     */
    void originalImageChanged(const QImage &image);

    /**
     * @brief Emitted when image matrix changes / Ispuskayetsya pri izmenenii matritsy izobrazheniya
     */
    void imageMatrixChanged();

    /**
     * @brief Emitted when attribute matrix changes / Ispuskayetsya pri izmenenii atributnoy matritsy
     */
    void attributeMatrixChanged();

    /**
     * @brief Emitted when selected edge changes / Ispuskayetsya pri izmenenii vybora granits
     */
    void selectedEdgeChanged(const QVector<QPoint>& edge);

    /**
     * @brief Emitted when sigma parameter changes / Ispuskayetsya pri izmenenii parametra sigma
     */
    void sigmaChanged(double sigma);

    /**
     * @brief Emitted when radius parameter changes / Ispuskayetsya pri izmenenii parametra radiusa
     */
    void radiusChanged(int radius);

    /**
     * @brief Emitted when edge selection mode toggles / Ispuskayetsya pri pereklyuchenii rezhima vybora granits
     */
    void edgeSelectionModeChanged(bool enabled);

    /**
     * @brief Emitted when any parameter changes (for batch updates) / Ispuskayetsya pri izmenenii lyubogo parametra (dlya gruppovykh obnovleniy)
     */
    void parametersChanged();

private:
    // ========================================================================
    // Private Helper Methods / Vspomogatel'nyye metody
    // ========================================================================

    /**
     * @brief Initializes matrices with correct dimensions from current image
     * @brief Initsializiruyet matritsy s korrektnymi razmerami iz tekushchego izobrazheniya
     */
    void initializeMatrices();

    /**
     * @brief Validates matrix dimensions against current image / Proveryayet razmery matrits s tekushchim izobrazheniyem
     */
    bool validateMatrices() const;

private:
    // ========================================================================
    // Image Data / Dannyye izobrazheniya
    // ========================================================================
    QImage m_currentImage;           // Current working image / Tekushcheye rabocheye izobrazheniye
    QImage m_originalImage;          // Original unmodified image / Original'noye neizmenennoye izobrazheniye

    // ========================================================================
    // Matrix Data / Matrichnyye dannyye
    // ========================================================================
    Matrix2D<double> m_imageMatrix;  // Grayscale image matrix / Matritsa izobrazheniya v ottenkakh serogo
    Matrix2D<int> m_attributeMatrix; // Edge attribute matrix / Atributnaya matritsa granits
    Matrix2D<double> m_testMatrix;   // Test matrix for profiles / Testovaya matritsa dlya profiley

    // ========================================================================
    // Selection Data / Dannyye vybora
    // ========================================================================
    QVector<QPoint> m_selectedEdge;  // Currently selected edge points / Vybrannyye tochki granits
    bool m_edgeSelectionMode = false; // Edge selection mode flag / Flag rezhima vybora granits

    // ========================================================================
    // Algorithm Parameters / Parametry algoritmov
    // ========================================================================
    double m_sigma = 5.0;            // Gaussian sigma / Sigma Gaussa
    int m_radius = 30;               // Kernel radius / Radius yadra
    double m_redMax = 255.0;         // Maximum red value / Maksimal'noye znacheniye krasnogo
    double m_blueMax = 255.0;        // Maximum blue value / Maksimal'noye znacheniye sinego
};

#endif // APPDATAMODEL_H
