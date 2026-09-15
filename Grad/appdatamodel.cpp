#include "appdatamodel.h"
#include "imageprocessor.h"
#include <QDebug>

// ============================================================================
// Constructor / Destructor / Konstruktor / Destruktor
// ============================================================================

/*
 * EN: Constructor - initializes empty data model
 * RU: Konstruktor - initsializiruyet pustuyu model' dannykh
 */
AppDataModel::AppDataModel(QObject *parent)
    : QObject(parent)
{
    // Initialize with empty data / Initsializatsiya s pustymi dannymi
}

/*
 * EN: Destructor - matrices are automatically cleaned up by std::vector
 * RU: Destruktor - matritsy avtomaticheski ochishchayutsya std::vector
 */
AppDataModel::~AppDataModel()
{
    // No manual cleanup needed / Ruchnaya ochistka ne trebuyetsya
    clear();
}

// ============================================================================
// Image Data Methods / Metody dannykh izobrazheniya
// ============================================================================

/*
 * EN: Sets current image and updates associated matrices
 * RU: Ustanavlivayet tekushcheye izobrazheniye i obnovlyayet svyazannyye matritsy
 */
void AppDataModel::setCurrentImage(const QImage &newCurrentImage)
{
    if (m_currentImage == newCurrentImage && !m_currentImage.isNull()) {
        return; // No change / Bez izmeneniy
    }

    m_currentImage = newCurrentImage;

    // Update matrices from new image / Obnovlyayem matritsy iz novogo izobrazheniya
    initializeMatrices();

    emit currentImageChanged(m_currentImage);
    emit imageMatrixChanged();
    emit attributeMatrixChanged();
}

/*
 * EN: Sets original image
 * RU: Ustanavlivayet original'noye izobrazheniye
 */
void AppDataModel::setOriginalImage(const QImage &newOriginalImage)
{
    if (m_originalImage == newOriginalImage && !m_originalImage.isNull()) {
        return;
    }

    m_originalImage = newOriginalImage;
    emit originalImageChanged(m_originalImage);
}

// ============================================================================
// Matrix Data Methods / Metody matrichnykh dannykh
// ============================================================================

/*
 * EN: Sets image matrix with validation
 * RU: Ustanavlivayet matritsu izobrazheniya s proverkoy
 */
void AppDataModel::setImageMatrix(const Matrix2D<double>& newImageMatrix)
{
    m_imageMatrix = newImageMatrix;
    emit imageMatrixChanged();
}

/*
 * EN: Sets attribute matrix with validation
 * RU: Ustanavlivayet atributnuyu matritsu s proverkoy
 */
void AppDataModel::setAttributeMatrix(const Matrix2D<int>& newAttributeMatrix)
{
    m_attributeMatrix = newAttributeMatrix;
    emit attributeMatrixChanged();
}

/*
 * EN: Sets test matrix
 * RU: Ustanavlivayet testovuyu matritsu
 */
void AppDataModel::setTestMatrix(const Matrix2D<double>& newTestMatrix)
{
    m_testMatrix = newTestMatrix;
}

// ============================================================================
// Edge Selection Methods / Metody vybora granits
// ============================================================================

/*
 * EN: Sets selected edge points
 * RU: Ustanavlivayet vybrannyye tochki granits
 */
void AppDataModel::setSelectedEdge(const QVector<QPoint>& newSelectedEdge)
{
    m_selectedEdge = newSelectedEdge;
    emit selectedEdgeChanged(m_selectedEdge);
}

/*
 * EN: Clears edge selection
 * RU: Ochishchayet vybor granits
 */
void AppDataModel::clearEdgeSelection()
{
    if (m_selectedEdge.isEmpty()) {
        return;
    }

    m_selectedEdge.clear();
    emit selectedEdgeChanged(m_selectedEdge);
}

/*
 * EN: Sets edge selection mode
 * RU: Ustanavlivayet rezhim vybora granits
 */
void AppDataModel::setEdgeSelectionMode(bool enabled)
{
    if (m_edgeSelectionMode == enabled) {
        return;
    }

    m_edgeSelectionMode = enabled;
    emit edgeSelectionModeChanged(m_edgeSelectionMode);

    if (!enabled) {
        // Clear selection when mode is disabled / Ochishchayem vybor pri otklyuchenii rezhima
        clearEdgeSelection();
    }
}

// ============================================================================
// Parameter Methods / Metody parametrov
// ============================================================================

/*
 * EN: Sets sigma value with bounds checking
 * RU: Ustanavlivayet znacheniye sigma s proverkoy granits
 */
void AppDataModel::setSigma(double newSigma)
{
    // First clamp to reasonable range / Snachala ogranichivayem v razumnykh predelakh
    newSigma = std::max(0.1, std::min(50.0, newSigma));

    if (std::abs(m_sigma - newSigma) < 1e-6) {
        return;
    }

    m_sigma = newSigma;
    emit sigmaChanged(m_sigma);
    emit parametersChanged();
}

/*
 * EN: Sets radius with bounds checking
 * RU: Ustanavlivayet radius s proverkoy granits
 */
void AppDataModel::setRadius(int newRadius)
{
    // Ensure radius is odd and positive / Obezpechivayem, chto radius nechyotnyy i polozhitel'nyy
    newRadius = std::max(1, newRadius);
    if (newRadius % 2 == 0) {
        newRadius += 1; // Make odd / Delayem nechyotnym
    }

    if (m_radius == newRadius) {
        return;
    }

    m_radius = newRadius;
    emit radiusChanged(m_radius);
    emit parametersChanged();
}

/*
 * EN: Sets maximum red value for visualization
 * RU: Ustanavlivayet maksimal'noye znacheniye krasnogo dlya vizualizatsii
 */
void AppDataModel::setRedMax(double newRedMax)
{
    newRedMax = std::max(0.0, std::min(255.0, newRedMax));

    if (std::abs(m_redMax - newRedMax) < 1e-6) {
        return;
    }

    m_redMax = newRedMax;
    emit parametersChanged();
}

/*
 * EN: Sets maximum blue value for visualization
 * RU: Ustanavlivayet maksimal'noye znacheniye sinego dlya vizualizatsii
 */
void AppDataModel::setBlueMax(double newBlueMax)
{
    newBlueMax = std::max(0.0, std::min(255.0, newBlueMax));

    if (std::abs(m_blueMax - newBlueMax) < 1e-6) {
        return;
    }

    m_blueMax = newBlueMax;
    emit parametersChanged();
}

// ============================================================================
// Utility Methods / Vspomogatel'nyye metody
// ============================================================================

/*
 * EN: Updates all matrices from current image
 * RU: Obnovlyayet vse matritsy iz tekushchego izobrazheniya
 */
void AppDataModel::updateMatricesFromCurrentImage()
{
    if (m_currentImage.isNull()) {
        qDebug() << "Warning: Cannot update matrices from null image";
        return;
    }

    // Convert current image to grayscale matrix / Preobrazuyem tekushcheye izobrazheniye v matritsu
    m_imageMatrix = ImageProcessor::fromGrayImage(m_currentImage);

    // Initialize attribute matrix with zeros / Initsializiruyem atributnuyu matritsu nulyami
    int width = m_currentImage.width();
    int height = m_currentImage.height();
    m_attributeMatrix = Matrix2D<int>(width, std::vector<int>(height, 0));

    // Initialize test matrix with zeros / Initsializiruyem testovuyu matritsu nulyami
    m_testMatrix = Matrix2D<double>(width, std::vector<double>(height, 0.0));

    emit imageMatrixChanged();
    emit attributeMatrixChanged();
}

/*
 * EN: Initializes matrices with correct dimensions from current image
 * RU: Initsializiruyet matritsy s korrektnymi razmerami iz tekushchego izobrazheniya
 */
void AppDataModel::initializeMatrices()
{
    if (m_currentImage.isNull()) {
        clear();
        return;
    }

    updateMatricesFromCurrentImage();
}

/*
 * EN: Validates that matrices have correct dimensions
 * RU: Proveryayet, chto matritsy imeyut korrektnyye razmery
 */
bool AppDataModel::validateMatrices() const
{
    if (m_currentImage.isNull()) {
        return false;
    }

    int width = m_currentImage.width();
    int height = m_currentImage.height();

    // Check image matrix / Proveryaem matritsu izobrazheniya
    if (static_cast<int>(m_imageMatrix.size()) != width ||
        (m_imageMatrix.size() > 0 && static_cast<int>(m_imageMatrix[0].size()) != height)) {
        qDebug() << "Warning: Image matrix dimension mismatch";
        return false;
    }

    // Check attribute matrix / Proveryaem atributnuyu matritsu
    if (static_cast<int>(m_attributeMatrix.size()) != width ||
        (m_attributeMatrix.size() > 0 && static_cast<int>(m_attributeMatrix[0].size()) != height)) {
        qDebug() << "Warning: Attribute matrix dimension mismatch";
        return false;
    }

    // Check test matrix / Proveryaem testovuyu matritsu
    if (static_cast<int>(m_testMatrix.size()) != width ||
        (m_testMatrix.size() > 0 && static_cast<int>(m_testMatrix[0].size()) != height)) {
        qDebug() << "Warning: Test matrix dimension mismatch";
        return false;
    }

    return true;
}

/*
 * EN: Clears all data
 * RU: Ochishchayet vse dannyye
 */
void AppDataModel::clear()
{
    m_currentImage = QImage();
    m_originalImage = QImage();

    // Clear matrices / Ochishchayem matritsy
    m_imageMatrix.clear();
    m_attributeMatrix.clear();
    m_testMatrix.clear();

    // Clear selection / Ochishchayem vybor
    m_selectedEdge.clear();
    m_edgeSelectionMode = false;

    // Emit signals / Ispuskayem signaly
    emit currentImageChanged(m_currentImage);
    emit imageMatrixChanged();
    emit attributeMatrixChanged();
    emit selectedEdgeChanged(m_selectedEdge);
    emit edgeSelectionModeChanged(false);
}

/*
 * EN: Resets current image to original
 * RU: Vosstanavlivayet tekushcheye izobrazheniye do original'nogo
 */
void AppDataModel::resetToOriginal()
{
    if (m_originalImage.isNull()) {
        qDebug() << "Warning: No original image to reset to";
        return;
    }

    setCurrentImage(m_originalImage);
}
