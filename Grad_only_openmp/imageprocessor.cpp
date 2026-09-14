#include "imageprocessor.h"
#include <QElapsedTimer>
#include <omp.h>

// ============================================================================
// Utility Functions / Vspomogatel'nye funkcii
// ============================================================================

/*
 * EN: Checks if vector contains point p
 * RU: Proveryayet, soderzhit li vektor tochku p
 */
bool ImageProcessor::containsVP(const QVector<QPoint>& vec, const QPoint& p)
{
    for(int i = 0; i < vec.size(); i++)
    {
        if(vec[i].x() == p.x() && vec[i].y() == p.y())
            return true;
    }
    return false;
}

/*
 * EN: If pos is an edge pixel, selects all connected edge pixels
 *     and updates the attribute matrix
 * RU: Yesli pos - eto kraevoy piksel', vybiraet vse svyazannye kraevye pikseli
 *     i obnovlyayet atributnuyu matricu
 */
void ImageProcessor::selectEdge(const QPoint& pos, Matrix2D<int>& iiImAttMat, QVector<QPoint>& selectedEdge)
{
    // Clear previous selection / Ochishchaem predydushchiy vybor
    selectedEdge.clear();

    // Check if starting point is an edge / Proveryaem, yavlyaetsya li nachal'naya tochka kraem
    if(iiImAttMat[pos.x()][pos.y()] != static_cast<int>(attribute::isEdge))
        return;

    selectedEdge.push_back(pos);
    iiImAttMat[pos.x()][pos.y()] = static_cast<int>(attribute::isSelectedEdge);

    size_t prevSize = 0, newSize = 1;

    // Flood fill to select all connected edge pixels / Zalivka dlya vybora vseh svyazannyh kraevyh pikseley
    while(prevSize < newSize)
    {
        prevSize = selectedEdge.size();
        for(int s = 0; s < selectedEdge.size(); s++)
        {
            const QPoint& current = selectedEdge.at(s);

            // Check 8-connected neighborhood / Proveryaem 8-svyaznuyu okrestnost'
            for(int ii = -1; ii <= 1; ii++)
            {
                for(int jj = -1; jj <= 1; jj++)
                {
                    // Skip the center point / Propuskaem tsentral'nuyu tochku
                    if(ii == 0 && jj == 0) continue;

                    int newX = current.x() + ii;
                    int newY = current.y() + jj;

                    // Bounds check / Proverka granits
                    if(newX < 0 || newX >= static_cast<int>(iiImAttMat.size()) ||
                       newY < 0 || newY >= static_cast<int>(iiImAttMat[0].size()))
                        continue;

                    if(iiImAttMat[newX][newY] == static_cast<int>(attribute::isEdge) &&
                       !containsVP(selectedEdge, QPoint(newX, newY)))
                    {
                        selectedEdge.push_back(QPoint(newX, newY));
                        iiImAttMat[newX][newY] = static_cast<int>(attribute::isSelectedEdge);
                    }
                }
            }
        }
        newSize = selectedEdge.size();
    }
}

// ============================================================================
// Profile Building Functions / Funkcii postroeniya profilya
// ============================================================================

/*
 * EN: Builds profile perpendicular to edge by sampling along gradient direction
 * RU: Stroit profil' perpendikulyarno krayu, vybiraya tochki vdol' gradienta
 */
ProfileResult ImageProcessor::buildProfile001(const ProfileParameters& params, const Matrix2D<double>& ddImMat)
{
    ProfileResult result;

    // Validate input / Proveryaem vhodnye dannye
    if(ddImMat.empty() || ddImMat[0].empty()) {
        qDebug() << "Error: Empty image matrix in buildProfile001";
        return result;
    }

    // Calculate gradient magnitude and normalize / Vychislyaem velichinu gradienta i normalizuem
    double gradientMagnitude = std::sqrt(params.xGradient * params.xGradient +
                                         params.yGradient * params.yGradient);

    if(gradientMagnitude < 1e-10) {
        qDebug() << "Warning: Zero gradient magnitude in buildProfile001";
        return result;
    }

    double normalizedXGrad = params.xGradient / gradientMagnitude;
    double normalizedYGrad = params.yGradient / gradientMagnitude;

    int numSteps = static_cast<int>(params.numSigma * params.sigma);
    if(numSteps <= 0) numSteps = 50;  // Default value / Znachenie po-umolchaniyu

    result.points.reserve(numSteps);
    result.indices.reserve(numSteps);

    int rows = static_cast<int>(ddImMat.size());
    int cols = static_cast<int>(ddImMat[0].size());

    for(int step = 0; step < numSteps; ++step) {
        // Calculate mu parameter / Vychislyaem parametr mu
        double muS = -params.numSigma * params.sigma +
                     2.0 * params.numSigma * params.sigma * step / numSteps;

        double x = params.x0 + muS * normalizedXGrad;
        double y = params.y0 + muS * normalizedYGrad;

        int xCoord = static_cast<int>(std::round(x));
        int yCoord = static_cast<int>(std::round(y));

        // Bounds check with clamping / Proverka granic s privedeniem
        if(xCoord < 0 || xCoord >= rows || yCoord < 0 || yCoord >= cols) {
            result.points.push_back(QPointF(static_cast<double>(step), 0.0));
            result.indices.push_back(QPoint(-1, -1));  // Invalid index / Nekorrektnyy indeks
        } else {
            double yProf = ddImMat[xCoord][yCoord];
            result.points.push_back(QPointF(static_cast<double>(step), yProf));
            result.indices.push_back(QPoint(xCoord, yCoord));
        }
    }
    return result;
}

/*
 * EN: Builds profile with Laplacian convolution for better edge detection
 * RU: Stroit profil' s primeneniyem svertki Laplasa dlya uluchshennogo vyyavleniya granic
 */
ProfileResult ImageProcessor::buildProfile002(const ProfileParameters& params,
                                                     const Matrix2D<double>& ddImMat,
                                                     int convRad)
{
    ProfileResult result;

    // Validate input / Proveryaem vhodnye dannye
    if(ddImMat.empty() || ddImMat[0].empty()) {
        qDebug() << "Error: Empty image matrix in buildProfile002";
        return result;
    }

    const int rows = static_cast<int>(ddImMat.size());
    const int cols = static_cast<int>(ddImMat[0].size());

    // Calculate gradient magnitude and normalize / Vychislyaem velichinu gradienta i normalizuem
    double gradientMagnitude = std::sqrt(params.xGradient * params.xGradient +
                                         params.yGradient * params.yGradient);

    if(gradientMagnitude < 1e-10) {
        qDebug() << "Error: Zero gradient magnitude in buildProfile002";
        return result;
    }

    double normalizedXGrad = params.xGradient / gradientMagnitude;
    double normalizedYGrad = params.yGradient / gradientMagnitude;

    int kernelSize = 2 * convRad + 1;
    Matrix2D<double> convCore = getLapl(kernelSize, kernelSize, params.sigma);

    if((int)convCore.size() != kernelSize || (int)convCore[0].size() != kernelSize) {
        qDebug() << "Error: Convolution kernel has wrong size";
        return result;
    }

    int numSteps = static_cast<int>(params.numSigma * params.sigma);
    if(numSteps <= 0) numSteps = 50;  // Default value / Znachenie po-umolchaniyu

    result.points.reserve(numSteps);
    result.indices.reserve(numSteps);

    // Pre-calculate kernel indices for faster access / Predvaritel'no vychislyaem indeksy yadra
    std::vector<int> kernelOffsets;
    kernelOffsets.reserve(kernelSize * kernelSize);
    for(int i = -convRad; i <= convRad; i++) {
        for(int j = -convRad; j <= convRad; j++) {
            kernelOffsets.push_back(i);
            kernelOffsets.push_back(j);
        }
    }

    for(int step = 0; step < numSteps; ++step) {
        double muS = -params.numSigma * params.sigma +
                     2.0 * params.numSigma * params.sigma * step / numSteps;

        double x = params.x0 + muS * normalizedXGrad;
        double y = params.y0 + muS * normalizedYGrad;

        int xCoord = static_cast<int>(std::round(x));
        int yCoord = static_cast<int>(std::round(y));

        if(xCoord < 0 || xCoord >= rows || yCoord < 0 || yCoord >= cols) {
            result.points.push_back(QPointF(static_cast<double>(step), 0.0));
            result.indices.push_back(QPoint(-1, -1));
        } else {
            double yProf = 0.0;

            // Apply Laplacian convolution / Primenyaem svertku Laplasa
            for(int ki = 0; ki < kernelSize; ki++) {
                int kernel_i = ki - convRad;
                for(int kj = 0; kj < kernelSize; kj++) {
                    int kernel_j = kj - convRad;

                    int convX = xCoord + kernel_i;
                    int convY = yCoord + kernel_j;

                    // Mirror boundary handling / Zerkal'naya obrabotka granits
                    if(convX < 0) convX = -convX;
                    else if(convX >= rows) convX = 2 * (rows - 1) - convX;

                    if(convY < 0) convY = -convY;
                    else if(convY >= cols) convY = 2 * (cols - 1) - convY;

                    // Final bounds check / Final'naya proverka granits
                    convX = std::max(0, std::min(convX, rows - 1));
                    convY = std::max(0, std::min(convY, cols - 1));

                    yProf += convCore[ki][kj] * ddImMat[convX][convY];
                }
            }

            result.points.push_back(QPointF(static_cast<double>(step), yProf));
            result.indices.push_back(QPoint(xCoord, yCoord));
        }
    }

    return result;
}

/*
 * EN: Cuts profile from center to nearest local min/max, then steps further
 * RU: Obrezayet profil' ot tsentra do blizhayshih lokal'nogo min i max, zatem otstupayet
 */
QVector<QPointF> ImageProcessor::cutProfile(const QVector<QPointF>& profile, int step)
{
    if(profile.size() < 3) return profile;

    QVector<QPointF> res;
    int left = profile.size()/2;
    int right = profile.size()/2;
    bool leftIncr = true;
    bool rightIncr = true;
    bool decrLeft = (profile[left-1].y() < profile[left].y());

    const int MAX_ITERATIONS = 1000;  // Prevent infinite loop / Predotvrashchaem beskonechnyy tsikl
    int iterations = 0;

    while((leftIncr || rightIncr) && iterations < MAX_ITERATIONS) {
        iterations++;

        if(left <= 0 || right >= profile.size()-1) break;

        if(decrLeft){
            if(profile[left-1].y() < profile[left].y() && leftIncr){
                left--;
            }
            else if(leftIncr){
                left = (left < step) ? 0 : left - step;
                leftIncr = false;
            }

            if(profile[right+1].y() > profile[right].y() && rightIncr){
                right++;
            }
            else if(rightIncr){
                right = (right > profile.size() - step - 1) ? profile.size() - 1 : right + step;
                rightIncr = false;
            }
        }
        else{
            if(profile[left-1].y() > profile[left].y() && leftIncr){
                left--;
            }
            else if(leftIncr){
                left = (left < step) ? 0 : left - step;
                leftIncr = false;
            }

            if(profile[right+1].y() < profile[right].y() && rightIncr){
                right++;
            }
            else if(rightIncr){
                right = (right > profile.size() - step - 1) ? profile.size() - 1 : right + step;
                rightIncr = false;
            }
        }
    }

    // Ensure valid range / Obezpechivaem korrektnyy diapazon
    left = std::max(0, left);
    right = std::min(static_cast<int>(profile.size() - 1), right);
    if(left >= right) return QVector<QPointF>();

    res.reserve(right - left);
    for(int i = left; i < right; i++){
        res.push_back(profile[i]);
    }
    return res;
}

/*
 * EN: Cuts profile between left and right indices with bounds checking
 * RU: Obrezayet profil' mezhdu levyym i pravym indeksami s proverkoy granits
 */
QVector<QPointF> ImageProcessor::cutProfile(const QVector<QPointF>& profile, int left, int right)
{
    QVector<QPointF> result;

    if(profile.empty()) return result;

    // Ensure valid bounds / Obezpechivaem korrektnye granitsy
    left = std::max(0, std::min(left, static_cast<int>(profile.size() - 1)));
    right = std::max(0, std::min(right, static_cast<int>(profile.size() - 1)));

    // Ensure left <= right / Ubezhdaemsya, chto left <= right
    if(left > right) std::swap(left, right);

    result.reserve(right - left + 1);
    for(int i = left; i <= right; i++){
        result.push_back(profile[i]);
    }
    return result;
}

/*
 * EN: Calculates RMS residual between two profiles
 * RU: Vychislyayet RMS nevyazku mezhdu dvumya profilyami
 */
double ImageProcessor::calculateResidualOfProfiles(const QVector<QPointF>& profile1,
                                                    const QVector<QPointF>& profile2)
{
    if(profile1.size() != profile2.size() || profile1.isEmpty()) {
        return std::numeric_limits<double>::max();
    }

    double sumSquared = 0.0;
    for(int i = 0; i < profile1.size(); ++i) {
        double diff = profile1[i].y() - profile2[i].y();
        sumSquared += diff * diff;
    }

    return std::sqrt(sumSquared / profile1.size());  // RMS, not just sum
}

/*
 * EN: Calculates sum of squared differences between profiles
 * RU: Vychislyayet summu kvadratov raznostey mezhdu profilyami
 */
double ImageProcessor::differenceOfTwoProfiles(const QVector<QPointF>& profile1,
                                                const QVector<QPointF>& profile2)
{
    if(profile1.size() != profile2.size()) return -1;

    double res = 0.0;
    for(int i = 0; i < profile1.size(); i++){
        double diff = profile1[i].y() - profile2[i].y();
        res += diff * diff;
    }
    return res;
}

/*
 * EN: Reinterpolates profile to new length using linear interpolation (full range)
 * RU: Pereinterpoliruyet profil' do novoy dliny s lineynoy interpolyatsiey (ves' diapazon)
 */
QVector<QPointF> ImageProcessor::reInterpolateProfile(const QVector<QPointF>& profile, int newLen)
{
    if(profile.empty() || newLen <= 0) return profile;
    return reInterpolateProfile(profile, 0.0, static_cast<double>(profile.size() - 1), newLen);
}

/*
 * EN: Reinterpolates profile from leftBound to rightBound to new length
 * RU: Pereinterpoliruyet profil' ot leftBound do rightBound do novoy dliny
 */
QVector<QPointF> ImageProcessor::reInterpolateProfile(const QVector<QPointF>& profile,
                                                       double leftBound,
                                                       double rightBound,
                                                       int newLen)
{
    QVector<QPointF> newProfile;

    // Validate inputs / Proveryaem vhodnye dannye
    if(profile.size() < 2 || newLen <= 0) {
        qDebug() << "Error: Invalid input in reInterpolateProfile";
        return profile;
    }

    // Clamp bounds to valid range / Ogranichivaem granitsy korrektnym diapazonom
    leftBound = std::max(0.0, std::min(leftBound, static_cast<double>(profile.size() - 1)));
    rightBound = std::max(0.0, std::min(rightBound, static_cast<double>(profile.size() - 1)));

    if(leftBound >= rightBound) {
        qDebug() << "Error: Invalid bounds in reInterpolateProfile";
        return profile;
    }

    // Special case: newLen == 1 / Osobyy sluchay: newLen == 1
    if(newLen == 1) {
        int midIdx = static_cast<int>((leftBound + rightBound) / 2.0);
        newProfile.append(QPointF(0.0, profile[midIdx].y()));
        return newProfile;
    }

    newProfile.reserve(newLen);
    double step = (rightBound - leftBound) / (newLen - 1);

    for(int i = 0; i < newLen; i++) {
        double currentPos = leftBound + i * step;

        int idx1 = static_cast<int>(std::floor(currentPos));
        int idx2 = static_cast<int>(std::ceil(currentPos));

        // Clamp indices to valid range / Ogranichivaem indeksy korrektnym diapazonom
        idx1 = std::max(0, std::min(idx1, static_cast<int>(profile.size() - 1)));
        idx2 = std::max(0, std::min(idx2, static_cast<int>(profile.size() - 1)));

        if(idx1 == idx2) {
            newProfile.append(QPointF(static_cast<double>(i), profile[idx1].y()));
            continue;
        }

        const QPointF& p1 = profile[idx1];
        const QPointF& p2 = profile[idx2];

        double t = (currentPos - idx1) / (idx2 - idx1);
        double y = p1.y() * (1.0 - t) + p2.y() * t;

        newProfile.append(QPointF(static_cast<double>(i), y));
    }

    return newProfile;
}

// ============================================================================
// Image Conversion Functions / Funkcii preobrazovaniya izobrazheniy
// ============================================================================

/*
 * EN: Converts QImage to grayscale double matrix
 * RU: Preobrazuyet QImage v seruyu matricu double
 */
Matrix2D<double> ImageProcessor::fromGrayImage(const QImage& image)
{
    if(image.isNull()) return Matrix2D<double>();

    int width = image.width();
    int height = image.height();
    Matrix2D<double> result(width, std::vector<double>(height));

    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            result[i][j] = static_cast<double>(qGray(image.pixel(i, j)));
        }
    }
    return result;
}

/*
 * EN: Converts double matrix to grayscale QImage with auto-scaling
 * RU: Preobrazuyet matricu double v seroye QImage s avtomaticheskim masshtabirovaniem
 */
QImage ImageProcessor::toGrayImage(const Matrix2D<double>& ddImageMat)
{
    if(ddImageMat.empty() || ddImageMat[0].empty()) {
        return QImage();
    }

    int width = static_cast<int>(ddImageMat.size());
    int height = static_cast<int>(ddImageMat[0].size());
    QImage resIm(width, height, QImage::Format_ARGB32);

    double dMin = min(ddImageMat);
    double dMax = max(ddImageMat);
    double range = dMax - dMin;

    // Avoid division by zero / Izbegayem deleniya na nol'
    if(range < 1e-10) range = 1.0;

    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            // Normalize to 0-255 range / Normalizuem v diapazon 0-255
            double normalized = (ddImageMat[i][j] - dMin) * 255.0 / range;
            int val = static_cast<int>(std::round(normalized));
            val = std::max(0, std::min(255, val));  // Clamp / Ogranichivaem
            resIm.setPixel(i, j, qRgb(val, val, val));
        }
    }
    return resIm;
}

/*
 * EN: Converts matrix to red/blue image (positive=red, negative=blue)
 * RU: Preobrazuyet matricu v krasno/sineye izobrazheniye (polozhitel'nye = krasnyy, otritsatel'nye = siniy)
 */
QImage ImageProcessor::toBlueRedImage(const Matrix2D<double>& ddImageMat, double dRedMax, double dBlueMax)
{
    if(ddImageMat.empty() || ddImageMat[0].empty()) {
        return QImage();
    }

    int width = static_cast<int>(ddImageMat.size());
    int height = static_cast<int>(ddImageMat[0].size());
    QImage resIm(width, height, QImage::Format_ARGB32);

    double dMin = min(ddImageMat);
    double dMax = max(ddImageMat);

    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            double val = ddImageMat[i][j];
            if(val < 0)
            {
                // Negative values become blue / Otritsatel'nye znacheniya stanovatsya sinimi
                int blueVal = static_cast<int>(std::round(-val * dBlueMax / std::max(-dMin, 1.0)));
                blueVal = std::max(0, std::min(255, blueVal));
                resIm.setPixel(i, j, qRgb(0, 0, blueVal));
            }
            else if(val > 0)
            {
                // Positive values become red / Polozhitel'nye znacheniya stanovatsya krasnymi
                int redVal = static_cast<int>(std::round(val * dRedMax / std::max(dMax, 1.0)));
                redVal = std::max(0, std::min(255, redVal));
                resIm.setPixel(i, j, qRgb(redVal, 0, 0));
            }
            else
            {
                resIm.setPixel(i, j, qRgb(0, 0, 0));
            }
        }
    }
    return resIm;
}

/*
 * EN: Combines image matrix with edge overlay (edges in white)
 * RU: Ob"yedinyayet matricu izobrazheniya s nalozheniyem granic (granicy belye)
 */
QImage ImageProcessor::combineImageWithEdge(const Matrix2D<double>& ddImageMat,
                                             const Matrix2D<int>& iiImageAttMat)
{
    if(ddImageMat.empty() || ddImageMat[0].empty()) return QImage();

    int width = static_cast<int>(ddImageMat.size());
    int height = static_cast<int>(ddImageMat[0].size());
    QImage resIm(width, height, QImage::Format_ARGB32);

    // Normalization for display / Normalizatsiya dlya otobrazheniya
    double dMin = min(ddImageMat);
    double dMax = max(ddImageMat);
    double range = dMax - dMin;
    if(range < 1e-10) range = 1.0;

    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            int attr = iiImageAttMat[i][j];

            if(attr == static_cast<int>(attribute::isEdge)) {
                resIm.setPixel(i, j, qRgb(255, 255, 255));  // White edge / Belaya granitsa
            }
            else if(attr == static_cast<int>(attribute::isSelectedEdge)) {
                resIm.setPixel(i, j, qRgb(0, 255, 255));    // Cyan selected / Golubaya vybrannaya
            }
            else {
                // Normalized grayscale / Normalizovannaya seraya shkala
                int val = static_cast<int>((ddImageMat[i][j] - dMin) * 255.0 / range);
                val = std::max(0, std::min(255, val));
                resIm.setPixel(i, j, qRgb(val, val, val));
            }
        }
    }
    return resIm;
}

/*
 * EN: Combines QImage with edge overlay
 * RU: Ob"yedinyayet QImage s nalozheniyem granic
 */
QImage ImageProcessor::combineImageWithEdge(const QImage& Image, const Matrix2D<int>& iiImageAttMat)
{
    if(Image.isNull()) return QImage();

    int width = Image.width();
    int height = Image.height();
    QImage resIm(width, height, QImage::Format_ARGB32);

    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            int attr = iiImageAttMat[i][j];

            if(attr == static_cast<int>(attribute::isEdge)) {
                resIm.setPixel(i, j, qRgb(255, 255, 255));
            }
            else if(attr == static_cast<int>(attribute::isSelectedEdge)) {
                resIm.setPixel(i, j, qRgb(0, 255, 255));
            }
            else {
                resIm.setPixel(i, j, Image.pixel(i, j));
            }
        }
    }
    return resIm;
}

/*
 * EN: Combines image with profile indices (colors them cyan)
 * RU: Ob"yedinyayet izobrazheniye s indeksami profilya (okrashivayet ikh v goluboy)
 */
QImage ImageProcessor::combineImageWithProfile(const QImage& im1, const QVector<QPoint>& profileIndices)
{
    if(im1.isNull()) return QImage();

    QImage resIm = im1;  // Start with copy / Nachinayem s kopii

    for(int i = 0; i < profileIndices.size(); i++){
        int x = profileIndices[i].x();
        int y = profileIndices[i].y();
        if(x >= 0 && x < im1.width() && y >= 0 && y < im1.height()) {
            resIm.setPixel(x, y, qRgb(0, 255, 255));
        }
    }
    return resIm;
}

// ============================================================================
// Convolution and Kernel Functions / Funkcii svertki i yader
// ============================================================================

/*
 * EN: Convolves QImage with kernel (RGB channels separately)
 * RU: Svertka QImage s yadrom (RGB kanaly otdel'no)
 */
QImage ImageProcessor::convImage(const QImage& image, const Matrix2D<double>& ddConvCore)
{
    if(image.isNull() || ddConvCore.empty()) return image;

    QImage resIm = image;
    int iRad = static_cast<int>(ddConvCore.size());
    int kernelRadius = iRad / 2;

    // Calculate kernel sum for normalization / Vychislyayem summu yadra dlya normalizatsii
    double kernelSum = 0.0;
    for(int ii = 0; ii < iRad; ii++) {
        for(int jj = 0; jj < iRad; jj++) {
            kernelSum += ddConvCore[ii][jj];
        }
    }
    if(std::abs(kernelSum) < 1e-10) kernelSum = 1.0;
    int width = image.width();
    int height = image.height();

// #pragma omp parallel for collapse(2)
    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            double rSum = 0, gSum = 0, bSum = 0;

            for(int ii = 0; ii < iRad; ii++)
            {
                int kernel_i = ii - kernelRadius;
                int it = i + kernel_i;

                // Mirror boundary handling / Zerkal'naya obrabotka granits
                if(it < 0) it = -it;
                else if(it >= width) it = 2 * (width - 1) - it;
                it = std::max(0, std::min(it, width - 1));

                for(int jj = 0; jj < iRad; jj++)
                {
                    int kernel_j = jj - kernelRadius;
                    int jt = j + kernel_j;

                    if(jt < 0) jt = -jt;
                    else if(jt >= height) jt = 2 * (height - 1) - jt;
                    jt = std::max(0, std::min(jt, height - 1));

                    double kernelVal = ddConvCore[ii][jj];
                    QRgb pixel = image.pixel(it, jt);

                    rSum += qRed(pixel) * kernelVal;
                    gSum += qGreen(pixel) * kernelVal;
                    bSum += qBlue(pixel) * kernelVal;
                }
            }

            rSum = std::max(0.0, std::min(255.0, rSum / kernelSum));
            gSum = std::max(0.0, std::min(255.0, gSum / kernelSum));
            bSum = std::max(0.0, std::min(255.0, bSum / kernelSum));

            resIm.setPixel(i, j, qRgb(static_cast<int>(rSum),
                                      static_cast<int>(gSum),
                                      static_cast<int>(bSum)));
        }
    }
    return resIm;
}

/*
 * EN: Convolves double matrix with kernel (optimized with precomputed indices)
 * RU: Svertka matricy double s yadrom (optimizirovano s predvychislennymi indeksami)
 */
Matrix2D<double> ImageProcessor::convMat(const Matrix2D<double>& ddImageMat,
                                          const Matrix2D<double>& ddConvCore)
{
    if(ddImageMat.empty() || ddImageMat[0].empty() || ddConvCore.empty()) {
        return Matrix2D<double>();
    }

    int rows = static_cast<int>(ddImageMat.size());
    int cols = static_cast<int>(ddImageMat[0].size());
    int kernelSize = static_cast<int>(ddConvCore.size());
    int kernelRadius = kernelSize / 2;

    Matrix2D<double> ddResMat(rows, std::vector<double>(cols, 0.0));

    // Calculate kernel sum for normalization / Vychislyayem summu yadra dlya normalizatsii
    double kernelSum = 0.0;
    for(int i = 0; i < kernelSize; i++) {
        for(int j = 0; j < kernelSize; j++) {
            kernelSum += ddConvCore[i][j];
        }
    }
    if(std::abs(kernelSum) < 1e-10) kernelSum = 1.0;

    // Precompute mirrored indices for faster access / Predvychislyaem zerkal'nye indeksy dlya bystrogo dostupa
    std::vector<int> xIndices(rows + kernelSize);
    std::vector<int> yIndices(cols + kernelSize);

    for(int i = -kernelRadius; i < rows + kernelRadius; i++) {
        int idx = i;
        if(i < 0) idx = -i;
        else if(i >= rows) idx = 2 * (rows - 1) - i;
        idx = std::max(0, std::min(idx, rows - 1));
        xIndices[i + kernelRadius] = idx;
    }

    for(int j = -kernelRadius; j < cols + kernelRadius; j++) {
        int idx = j;
        if(j < 0) idx = -j;
        else if(j >= cols) idx = 2 * (cols - 1) - j;
        idx = std::max(0, std::min(idx, cols - 1));
        yIndices[j + kernelRadius] = idx;
    }

    // Perform convolution / Vypolnyayem svertku
#pragma omp parallel for collapse(2)
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            double matSum = 0.0;

            for(int ki = 0; ki < kernelSize; ki++) {
                int imgX = xIndices[i + ki];

                for(int kj = 0; kj < kernelSize; kj++) {
                    int imgY = yIndices[j + kj];
                    matSum += ddImageMat[imgX][imgY] * ddConvCore[ki][kj];
                }
            }

            ddResMat[i][j] = matSum / kernelSum;
        }
    }

    return ddResMat;
}

/*
 * EN: Generates Gaussian kernel (normalized)
 * RU: Generiruyet yadro Gaussa (normalizovannoye)
 */
Matrix2D<double> ImageProcessor::getGauss(int iXsize, int iYsize, double dSigma)
{
    if(iXsize <= 0 || iYsize <= 0 || dSigma <= 0) {
        return Matrix2D<double>();
    }

    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize, 0.0));
    int xCenter = iXsize / 2;
    int yCenter = iYsize / 2;
    double sigmaSq = dSigma * dSigma;
    double twoSigmaSq = 2.0 * sigmaSq;
    double sum = 0.0;

    for(int i = 0; i < iXsize; i++)
    {
        int dx = i - xCenter;
        int dxSq = dx * dx;

        for(int j = 0; j < iYsize; j++)
        {
            int dy = j - yCenter;
            double radSq = dxSq + dy * dy;
            ddRes[i][j] = exp(-radSq / twoSigmaSq);
            sum += ddRes[i][j];
        }
    }

    // Normalize so sum = 1 / Normalizuyem tak, chtoby summa = 1
    if(sum > 0) {
        for(int i = 0; i < iXsize; i++) {
            for(int j = 0; j < iYsize; j++) {
                // ddRes[i][j] /= sum;
            }
        }
    }

    qDebug() << "Gauss kernel sum =" << sum;
    return ddRes;
}

/*
 * EN: Generates X-gradient kernel (Gaussian derivative)
 * RU: Generiruyet yadro X-gradienta (proizvodnaya Gaussa)
 */
Matrix2D<double> ImageProcessor::getXGradCore(int iXsize, int iYsize, double dSigma)
{
    if(iXsize <= 0 || iYsize <= 0 || dSigma <= 0) {
        return Matrix2D<double>();
    }

    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize, 0.0));
    int xCenter = iXsize / 2;
    int yCenter = iYsize / 2;
    double sigmaSq = dSigma * dSigma;
    double twoSigmaSq = 2.0 * sigmaSq;
    double sum = 0.0;

    for(int i = 0; i < iXsize; i++)
    {
        int dx = i - xCenter;
        int dxSq = dx * dx;

        for(int j = 0; j < iYsize; j++)
        {
            int dy = j - yCenter;
            double radSq = dxSq + dy * dy;
            ddRes[i][j] = static_cast<double>(dx) / sigmaSq * exp(-radSq / twoSigmaSq);
            sum += ddRes[i][j];
        }
    }

    qDebug() << "X grad kernel sum =" << sum;
    return ddRes;
}

/*
 * EN: Generates Y-gradient kernel (Gaussian derivative)
 * RU: Generiruyet yadro Y-gradienta (proizvodnaya Gaussa)
 */
Matrix2D<double> ImageProcessor::getYGradCore(int iXsize, int iYsize, double dSigma)
{
    if(iXsize <= 0 || iYsize <= 0 || dSigma <= 0) {
        return Matrix2D<double>();
    }

    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize, 0.0));
    int xCenter = iXsize / 2;
    int yCenter = iYsize / 2;
    double sigmaSq = dSigma * dSigma;
    double twoSigmaSq = 2.0 * sigmaSq;
    double sum = 0.0;

    for(int i = 0; i < iXsize; i++)
    {
        int dx = i - xCenter;
        int dxSq = dx * dx;

        for(int j = 0; j < iYsize; j++)
        {
            int dy = j - yCenter;
            double radSq = dxSq + dy * dy;
            ddRes[i][j] = static_cast<double>(dy) / sigmaSq * exp(-radSq / twoSigmaSq);
            sum += ddRes[i][j];
        }
    }

    qDebug() << "Y grad kernel sum =" << sum;
    return ddRes;
}

/*
 * EN: Generates Laplacian of Gaussian (LoG) kernel
 * RU: Generiruyet yadro Laplasiana Gaussa (LoG)
 */
Matrix2D<double> ImageProcessor::getLapl(int iXsize, int iYsize, double dSigma)
{
    if(iXsize <= 0 || iYsize <= 0 || dSigma <= 0) {
        return Matrix2D<double>();
    }

    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize, 0.0));
    int xCenter = iXsize / 2;
    int yCenter = iYsize / 2;
    double sigmaSq = dSigma * dSigma;
    double twoSigmaSq = 2.0 * sigmaSq;
    double sum = 0.0;

    for(int i = 0; i < iXsize; i++)
    {
        int dx = i - xCenter;
        int dxSq = dx * dx;

        for(int j = 0; j < iYsize; j++)
        {
            int dy = j - yCenter;
            double radSq = dxSq + dy * dy;
            // LoG: (r^2/sigma^2 - 2) * exp(-r^2/(2*sigma^2))
            ddRes[i][j] = (radSq / sigmaSq - 2.0) * exp(-radSq / twoSigmaSq);
            sum += ddRes[i][j];
        }
    }

    qDebug() << "Laplacian kernel sum =" << sum;
    return ddRes;
}

// ============================================================================
// Edge Detection Functions / Funkcii detekcii granic
// ============================================================================

/*
 * EN: Finds edge pixels where positive and negative neighbors exist in radius iRad
 * RU: Nahodit kraevye pikseli, gde sushchestvuyut polozhitel'nye i otritsatel'nye sosedi v radiuse iRad
 */
Matrix2D<double> ImageProcessor::findEdges(const Matrix2D<double>& ddImageMat,
                                            int iRad,
                                            Matrix2D<int>& iiImageAttMat)
{
    if(ddImageMat.empty() || ddImageMat[0].empty()) return Matrix2D<double>();

    int rows = static_cast<int>(ddImageMat.size());
    int cols = static_cast<int>(ddImageMat[0].size());
    Matrix2D<double> ddImageEdgeMat(rows, std::vector<double>(cols, 0.0));

    int kernelRadius = iRad / 2;

#pragma omp parallel for collapse(2)
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            bool bPos = false, bNeg = false;

            for(int ii = 0; ii < iRad; ii++)
            {
                int kernel_i = ii - kernelRadius;
                int it = i + kernel_i;

                // Mirror boundary / Zerkal'naya granitsa
                if(it < 0) it = -it;
                else if(it >= rows) it = 2 * (rows - 1) - it;
                it = std::max(0, std::min(it, rows - 1));

                for(int jj = 0; jj < iRad; jj++)
                {
                    int kernel_j = jj - kernelRadius;
                    int jt = j + kernel_j;

                    if(jt < 0) jt = -jt;
                    else if(jt >= cols) jt = 2 * (cols - 1) - jt;
                    jt = std::max(0, std::min(jt, cols - 1));

                    double val = ddImageMat[it][jt];
                    if(val > 1e-10) bPos = true;
                    else if(val < -1e-10) bNeg = true;

                    if(bPos && bNeg) break;
                }
                if(bPos && bNeg) break;
            }

            if(bPos && bNeg)
            {
                iiImageAttMat[i][j] = static_cast<int>(attribute::isEdge);
                ddImageEdgeMat[i][j] = 255.0;
            }
            else
            {
                iiImageAttMat[i][j] = 0;
                ddImageEdgeMat[i][j] = 0.0;
            }
        }
    }
    return ddImageEdgeMat;
}

/*
 * EN: Finds edge pixels without attribute matrix (simplified version)
 * RU: Nahodit kraevye pikseli bez atributnoy matritsy (upreshchennaya versiya)
 */
Matrix2D<double> ImageProcessor::findEdges(const Matrix2D<double>& ddImageMat, int iRad)
{
    if(ddImageMat.empty() || ddImageMat[0].empty()) return Matrix2D<double>();

    int rows = static_cast<int>(ddImageMat.size());
    int cols = static_cast<int>(ddImageMat[0].size());
    Matrix2D<double> ddImageEdgeMat(rows, std::vector<double>(cols, 0.0));

    int kernelRadius = iRad / 2;

#pragma omp parallel for collapse(2)
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            bool bPos = false, bNeg = false;

            for(int ii = 0; ii < iRad; ii++)
            {
                int kernel_i = ii - kernelRadius;
                int it = i + kernel_i;

                if(it < 0) it = -it;
                else if(it >= rows) it = 2 * (rows - 1) - it;
                it = std::max(0, std::min(it, rows - 1));

                for(int jj = 0; jj < iRad; jj++)
                {
                    int kernel_j = jj - kernelRadius;
                    int jt = j + kernel_j;

                    if(jt < 0) jt = -jt;
                    else if(jt >= cols) jt = 2 * (cols - 1) - jt;
                    jt = std::max(0, std::min(jt, cols - 1));

                    double val = ddImageMat[it][jt];
                    if(val > 0) bPos = true;
                    else if(val < 0) bNeg = true;

                    if(bPos && bNeg) break;
                }
                if(bPos && bNeg) break;
            }

            if(bPos && bNeg)
            {
                ddImageEdgeMat[i][j] = 255.0;
            }
        }
    }
    return ddImageEdgeMat;
}

/*
 * EN: Finds edges with red/blue encoding (positive=red, negative=blue)
 * RU: Nahodit granicy s krasno/sinim kodirovaniyem (polozhitel'nye = krasnyy, otritsatel'nye = siniy)
 */
Matrix2D<double> ImageProcessor::findEdgesRB(const Matrix2D<double>& ddImageMat,
                                              int iRad,
                                              Matrix2D<int>& iiImageAttMat)
{
    if(ddImageMat.empty() || ddImageMat[0].empty()) return Matrix2D<double>();

    int rows = static_cast<int>(ddImageMat.size());
    int cols = static_cast<int>(ddImageMat[0].size());
    Matrix2D<double> ddImageEdgeMat(rows, std::vector<double>(cols, 0.0));

    int kernelRadius = iRad / 2;

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            bool bPos = false, bNeg = false;

            for(int ii = 0; ii < iRad; ii++)
            {
                int kernel_i = ii - kernelRadius;
                int it = i + kernel_i;

                if(it < 0) it = -it;
                else if(it >= rows) it = 2 * (rows - 1) - it;
                it = std::max(0, std::min(it, rows - 1));

                for(int jj = 0; jj < iRad; jj++)
                {
                    int kernel_j = jj - kernelRadius;
                    int jt = j + kernel_j;

                    if(jt < 0) jt = -jt;
                    else if(jt >= cols) jt = 2 * (cols - 1) - jt;
                    jt = std::max(0, std::min(jt, cols - 1));

                    double val = ddImageMat[it][jt];
                    if(val > 0) bPos = true;
                    else if(val < 0) bNeg = true;

                    if(bPos && bNeg) break;
                }
                if(bPos && bNeg) break;
            }

            if(bPos && bNeg)
            {
                // Positive: +1, Negative: -1 / Polozhitel'nye: +1, Otritsatel'nye: -1
                ddImageEdgeMat[i][j] = (ddImageMat[i][j] > 0) ? 1.0 : -1.0;
                iiImageAttMat[i][j] = static_cast<int>(attribute::isEdge);
            }
            else
            {
                ddImageEdgeMat[i][j] = 0.0;
                iiImageAttMat[i][j] = 0;
            }
        }
    }
    return ddImageEdgeMat;
}

/*
 * EN: Gaussian edge detection with attribute matrix (difference of Gaussians)
 * RU: Detektsiya granic po Gaussu s atributnoy matritsey (raznost' Gausianov)
 */
QImage ImageProcessor::gaussianEdgeDetection(const Matrix2D<double>& ddImageMat,
                                              double dSigma,
                                              int iRad,
                                              Matrix2D<int>& iiImageAttMat)
{

    if(ddImageMat.empty() || ddImageMat[0].empty()) return QImage();

    // Difference of Gaussians (DoG) / Raznost' Gausianov
    Matrix2D<double> gauss1 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma));
    Matrix2D<double> gauss2 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma * 1.6));
    Matrix2D<double> dog = elementWiseOperation(gauss1, gauss2, MatrixLambdas::Subtract<double>{});

    // Find zero crossings / Nahodim perekhody cherez nol'
    Matrix2D<double> edges = findEdges(dog, 3, iiImageAttMat);

    return toGrayImage(edges);
}

/*
 * EN: Simplified Gaussian edge detection (no attribute matrix)
 * RU: Uproshchennaya detektsiya granic po Gaussu (bez atributnoy matritsy)
 */
QImage ImageProcessor::gaussianEdgeDetection(const Matrix2D<double>& ddImageMat,
                                              double dSigma,
                                              int iRad)
{
    if(ddImageMat.empty() || ddImageMat[0].empty()) return QImage();

    // Difference of Gaussians (DoG) with different sigma ratio / Raznost' Gausianov s drugim otnosheniem sigma
    Matrix2D<double> gauss1 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma * 2.0));
    Matrix2D<double> gauss2 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma));
    Matrix2D<double> dog = elementWiseOperation(gauss1, gauss2, MatrixLambdas::Subtract<double>{});

    Matrix2D<double> edges = findEdges(dog, 3);

    return toGrayImage(edges);
}

// ============================================================================
// Sample Image Generation / Generatsiya testovyh izobrazheniy
// ============================================================================

/*
 * EN: Creates test image with two hollow circles (small version 200x200)
 * RU: Sozdayot testovoye izobrazheniye s dvumya pustotami (malaya versiya 200x200)
 */
QImage ImageProcessor::sampleTwoHollows()
{
    QImage res(200, 200, QImage::Format_ARGB32);
    res.fill(Qt::black);

    // First circle (bottom)
    int centerX1 = 100, centerY1 = 135, radius1 = 40;
    // Second circle (top)
    int centerX2 = 100, centerY2 = 65, radius2 = 30;

    for(int i = 0; i < 200; i++){
        for(int j = 0; j < 200; j++){
            int val = 0;

            int dx1 = i - centerX1;
            int dy1 = j - centerY1;
            double dist1 = sqrt(dx1*dx1 + dy1*dy1);

            if(dist1 < radius1){
                double intensity = sqrt(radius1 * radius1 - dist1 * dist1);
                intensity = std::max(0.0, std::min(255.0, intensity));
                val = std::max(val, static_cast<int>(intensity));
            }

            int dx2 = i - centerX2;
            int dy2 = j - centerY2;
            double dist2 = sqrt(dx2*dx2 + dy2*dy2);

            if(dist2 < radius2){
                double intensity = sqrt(radius2 * radius2 - dist2 * dist2);
                intensity = std::max(0.0, std::min(255.0, intensity));
                val = std::max(val, static_cast<int>(intensity));
            }


            res.setPixel(i, j, qRgb(val, val, val));
        }
    }
    return res;
}

/*
 * EN: Creates test image with two hollow circles (large version 400x550)
 * RU: Sozdayot testovoye izobrazheniye s dvumya pustotami (bol'shaya versiya 400x550)
 */
QImage ImageProcessor::sampleTwoHollowsBig()
{
    QImage res(400, 550, QImage::Format_ARGB32);
    res.fill(Qt::black);

    // First circle (center) / Pervaya okruzhnost' (tsentr)
    int centerX1 = 200, centerY1 = 200, radius1 = 100;
    // Second circle (lower) / Vtoraya okruzhnost' (nizhe)
    int centerX2 = 200, centerY2 = 375, radius2 = 75;

    for(int i = 0; i < 400; i++){
        for(int j = 0; j < 550; j++){
            // Check first circle / Proveryaem pervuyu okruzhnost'
            int dx1 = i - centerX1;
            int dy1 = j - centerY1;
            int distSq1 = dx1*dx1 + dy1*dy1;

            if(distSq1 < radius1 * radius1){
                int val = static_cast<int>(std::sqrt(radius1*radius1 - distSq1));
                val = std::max(0, std::min(255, val));
                res.setPixel(i, j, qRgb(val, val, val));
            }

            // Check second circle / Proveryaem vtoruyu okruzhnost'
            int dx2 = i - centerX2;
            int dy2 = j - centerY2;
            int distSq2 = dx2*dx2 + dy2*dy2;

            if(distSq2 < radius2 * radius2){
                int val = static_cast<int>(std::sqrt(radius2*radius2 - distSq2));
                val = std::max(0, std::min(255, val));
                res.setPixel(i, j, qRgb(val, val, val));
            }
        }
    }
    return res;
}

// ============================================================================
// Edge Refinement Functions / Funktsii utochneniya granits
// ============================================================================

/*
 * EN: Refines a single edge point using PCR (Profile Correlation Refinement) algorithm
 * RU: Utochnyayet odnu granichnuyu tochku s ispol'zovaniyem algoritma PCR (Profile Correlation Refinement)
 *
 * EN: This is the main edge refinement algorithm that uses profile correlation
 *     to find sub-pixel accurate edge positions
 * RU: Eto osnovnoy algoritm utochneniya granits, ispol'zuyushchiy korrelyatsiyu profiley
 *     dlya poiska subpiksel'nykh pozitsiy granits
 */
RefinementResult ImageProcessor::refineSinglePoint001(int n0, int m0, const RefinementParameters& params)
{
    RefinementResult result;

    // Initialize with default values / Initsializatsiya so znacheniyami po-umolchaniyu
    result.success = false;
    result.FFF1_0 = 0.0;
    result.FFF1_1 = 0.0;
    result.residual = std::numeric_limits<double>::max();
    result.refinedPosition = QPointF(static_cast<double>(n0), static_cast<double>(m0));

    // Bounds check / Proverka granits
    if(n0 < 0 || n0 >= params.NX || m0 < 0 || m0 >= params.NY) {
        qDebug() << "Error: Point out of bounds / Oshibka: Tochka vne granits";
        return result;
    }

    // Constants for profile building / Konstanty dlya postroeniya profilya
    const int n_sigma = params.n_sigma;
    const int n_myu = params.n_myu;
    const int x0 = n0, y0 = m0;
    const double sigma_myu = std::min(params.sigma1, params.sigma2);

    // Pre-allocate vectors / Predvaritel'no vydelyaem pamyat'
    QVector<double> mu;
    QVector<double> prof1, prof2;
    mu.reserve(n_myu);
    prof1.reserve(n_myu);
    prof2.reserve(n_myu);

    // Generate mu values / Generiruyem znacheniya mu
    double muStep = 2.0 * n_sigma * sigma_myu / n_myu;
    double muStart = -n_sigma * sigma_myu;
    for(int s = 0; s < n_myu; s++) {
        mu.push_back(muStart + s * muStep);
    }

    // Pre-calculate Laplacian kernel for efficiency / Predvaritel'no vychislyaem yadro Laplasa dlya effektivnosti
    // Note: This could be optimized further by using separable filters
    // Primechaniye: Mozhno dalee optimizirovat', ispol'zuya razdelyayemyye fil'try

    // Build profiles along gradient direction / Stroim profili vdol' napravleniya gradienta
    // Pre-allocate vectors for thread safety / Predvaritel'no vydelyaem vektory dlya potokobezopasnosti
    prof1.resize(n_myu);
    prof2.resize(n_myu);

// #pragma omp parallel for
    for(int s = 0; s < n_myu; s++) {
        double prof1Val = 0.0;
        double prof2Val = 0.0;
        double mu_s = mu[s];

        for(int n = 0; n < params.NX; n++) {
            double dx_base = x0 - n;

            for(int m = 0; m < params.NY; m++) {
                double dy_base = y0 - m;

                double dx = dx_base + mu_s * params.ex;
                double dy = dy_base + mu_s * params.ey;
                double radSq = dx*dx + dy*dy;

                double mult = (radSq / (params.sigma1 * params.sigma1) - 2.0) *
                              exp(-0.5 * radSq / (params.sigma1 * params.sigma1));

                prof1Val += params.A[n][m] * mult;
                prof2Val += params.B01[n][m] * mult;
            }
        }
        prof1[s] = prof1Val;
        prof2[s] = prof2Val;
    }

    // Find zero crossings / Nakhodim perekhody cherez nol'
    App_Stats y2_stats;
    y2_stats.gather_stats(prof2);
    int nmumax = std::max(y2_stats.x_max, y2_stats.x_min);
    int nmumin = std::min(y2_stats.x_max, y2_stats.x_min);

    int nL_Zero = n_myu / 2;
    int nR_Zero = n_myu / 2;
    for(int n = nmumin; n < nmumax && n < n_myu - 1; n++) {
        if(prof1[n] * prof1[n+1] < 0) nL_Zero = n;
        if(prof2[n] * prof2[n+1] < 0) nR_Zero = n;
    }

    int N_Zero = (nL_Zero + nR_Zero) / 2;

    // Constants for profile bound detection / Konstanty dlya poiska granits profilya
    const int MAX_COUNTER = 5;      // Number of consecutive increases/decreases to detect plateau
    const int MIN_MARGIN = 20;      // Minimum margin from edges
    const int SAFE_MARGIN = 35;     // Safe margin to avoid boundary issues

    // Find profile bounds by counting increases/decreases
    // Nakhodim granitsy profilya, schitaya uvelicheniya/umen'sheniya
    int XL1 = 0, XR1 = 0;
    int decrCounter = 0, incrCounter = 0;

    // Search right bound / Poisk pravoy granitsy
    for(int i = N_Zero + 1; i < n_myu - 1; i++) {
        if(prof2[i] < prof2[i+1]) {
            incrCounter = std::min(incrCounter + 1, MAX_COUNTER);
        }
        if(prof2[i] > prof2[i+1]) {
            decrCounter = std::min(decrCounter + 1, MAX_COUNTER);
        }
        if(decrCounter == MAX_COUNTER && incrCounter == MAX_COUNTER) {
            XR1 = std::min(i, n_myu - SAFE_MARGIN);
            break;
        }
    }
    if(XR1 == 0) XR1 = std::min(N_Zero + 50, n_myu - SAFE_MARGIN);

    // Search left bound / Poisk levoy granitsy
    decrCounter = 0;
    incrCounter = 0;
    for(int i = N_Zero - 1; i > 1; i--) {
        if(prof2[i] < prof2[i+1]) {
            incrCounter = std::min(incrCounter + 1, MAX_COUNTER);
        }
        if(prof2[i] > prof2[i+1]) {
            decrCounter = std::min(decrCounter + 1, MAX_COUNTER);
        }
        if(decrCounter == MAX_COUNTER && incrCounter == MAX_COUNTER) {
            XL1 = std::max(i, MIN_MARGIN);
            break;
        }
    }
    if(XL1 == 0) XL1 = std::max(N_Zero - 50, MIN_MARGIN);

    // Prepare first profile for interpolation / Podgotavlivayem pervyy profil' dlya interpolyatsii
    QVector<QPointF> yP1;
    yP1.reserve(XR1 - XL1);
    for(int i = XL1; i < XR1; i++) {
        yP1.push_back({static_cast<double>(i), prof1[i]});
    }

    const int NN = params.NN;

    // Reinterpolate first profile / Pereinterpoliruyem pervyy profil'
    QVector<QPointF> yyP1 = reInterpolateProfile(yP1, NN);
    QVector<double> yy1;
    yy1.reserve(NN);
    for(int i = 0; i < NN; i++) {
        yy1.push_back(yyP1[i].y());
    }

    // Search for optimal second profile bounds / Poisk optimal'nykh granits vtorogo profilya
    const int otstup = params.otstup;
    double minraz = std::numeric_limits<double>::max();
    int best_XL2 = XL1, best_XR2 = XR1;

    // Prepare full second profile for interpolation / Podgotavlivayem polnyy vtoroy profil' dlya interpolyatsii
    QVector<QPointF> yP2_full;
    yP2_full.reserve(n_myu);
    for(int i = 0; i < n_myu; i++) {
        yP2_full.push_back({static_cast<double>(i), prof2[i]});
    }

    // Search over possible bounds / Poisk po vozmozhnym granitsam
// #pragma omp parallel for collapse(2)
    for(int XL2 = XL1; XL2 > XL1 - otstup; XL2--) {
        for(int XR2 = XR1; XR2 < XR1 + otstup; XR2++) {
            if(XL2 < 0 || XR2 >= n_myu || XL2 >= XR2) continue;

            QVector<QPointF> yyP2 = reInterpolateProfile(yP2_full, XL2, XR2, NN);
            QVector<double> yy2;
            yy2.reserve(NN);
            for(int i = 0; i < NN; i++) {
                yy2.push_back(yyP2[i].y());
            }

            // Normalize both profiles / Normalizuyem oba profilya
            double m1 = 0.0, m2 = 0.0;
            for(int i = 0; i < NN; i++) {
                m1 += yy1[i];
                m2 += yy2[i];
            }
            m1 /= NN;
            m2 /= NN;

            double D1 = 0.0, D2 = 0.0;
            for(int i = 0; i < NN; i++) {
                D1 += (yy1[i] - m1) * (yy1[i] - m1);
                D2 += (yy2[i] - m2) * (yy2[i] - m2);
            }
            D1 = std::sqrt(D1 / NN);
            D2 = std::sqrt(D2 / NN);

            if(D1 < 1e-10) D1 = 1.0;
            if(D2 < 1e-10) D2 = 1.0;

            double MINRAZ = 0.0;
            for(int i = 0; i < NN; i++) {
                double norm1 = (yy1[i] - m1) / D1;
                double norm2 = (yy2[i] - m2) / D2;
                double diff = norm1 - norm2;
                MINRAZ += diff * diff;
            }

// #pragma omp critical
            {
                if(minraz > MINRAZ) {
                    minraz = MINRAZ;
                    best_XL2 = XL2;
                    best_XR2 = XR2;
                }
            }
        }
    }

    // Final processing with optimal bounds / Final'naya obrabotka s optimal'nymi granitsami
    QVector<QPointF> yyP2_final = reInterpolateProfile(yP2_full, best_XL2, best_XR2, NN);
    QVector<double> yy2_final;
    yy2_final.reserve(NN);
    for(int i = 0; i < NN; i++) {
        yy2_final.push_back(yyP2_final[i].y());
    }

    // Final normalization / Final'naya normalizatsiya
    double m1_f = 0.0, m2_f = 0.0;
    for(int i = 0; i < NN; i++) {
        m1_f += yy1[i];
        m2_f += yy2_final[i];
    }
    m1_f /= NN;
    m2_f /= NN;

    double D1_f = 0.0, D2_f = 0.0;
    for(int i = 0; i < NN; i++) {
        D1_f += (yy1[i] - m1_f) * (yy1[i] - m1_f);
        D2_f += (yy2_final[i] - m2_f) * (yy2_final[i] - m2_f);
    }
    D1_f = std::sqrt(D1_f / NN);
    D2_f = std::sqrt(D2_f / NN);

    if(D1_f < 1e-10) D1_f = 1.0;
    if(D2_f < 1e-10) D2_f = 1.0;

    // Calculate final residual / Vychislyayem final'nuyu nevyazku
    double MINRAZ_final = 0.0;
    for(int i = 0; i < NN; i++) {
        double norm1 = (yy1[i] - m1_f) / D1_f;
        double norm2 = (yy2_final[i] - m2_f) / D2_f;
        double diff = norm1 - norm2;
        MINRAZ_final += diff * diff;
    }

    // Calculate shift and scale coefficients / Vychislyayem koeffitsienty sdviga i masshtaba
    double xa = static_cast<double>(XL1);
    double xb = static_cast<double>(XR1);
    double ya = static_cast<double>(best_XL2 - XL1);
    double yb = static_cast<double>(best_XR2 - XR1);

    // Avoid division by zero / Izbegayem deleniya na nol'
    if(std::abs(yb - ya) < 1e-10) {
        qDebug() << "Error: Division by zero in shift calculation";
        return result;
    }

    // Find intersection point / Nakhodim tochku peresecheniya
    double x00 = (xa * yb - xb * ya) / (yb - ya);

    // Convert to mu space / Preobrazuyem v prostranstvo mu
    double mu00 = -n_sigma * sigma_myu + x00 * (2.0 * n_sigma * sigma_myu) / n_myu;

    // Final shift and scale coefficients / Final'nyye koeffitsienty sdviga i masshtaba
    result.FFF1_1 = mu00;      // Shift coefficient / Koeffitsient sdviga
    result.FFF1_0 = (yb - ya) / (xb - xa);  // Scale coefficient / Koeffitsient masshtaba
    result.residual = MINRAZ_final;

    // Calculate refined position / Vychislyayem utochnennuyu pozitsiyu
    double n_new = n0 + params.ex * result.FFF1_1;
    double m_new = m0 + params.ey * result.FFF1_1;
    result.refinedPosition = QPointF(n_new, m_new);
    result.success = true;

    return result;
}

// Add these to imageprocessor.cpp after the existing functions

// ============================================================================
// Profile Building Between Points / Postroeniye profilya mezhdu tochkami
// ============================================================================

/*
 * EN: Builds profile along a line between two points using Bresenham's algorithm
 * RU: Stroit profil' vdol' linii mezhdu dvumya tochkami s ispol'zovaniyem algoritma Brezenhema
 */
QVector<double> ImageProcessor::buildProfileBetweenPoints(const QPoint& p1, const QPoint& p2,
                                                           const QImage& image)
{
    QVector<double> profile;

    if (image.isNull() || p1 == p2) {
        return profile;
    }

    int x1 = p1.x(), y1 = p1.y();
    int x2 = p2.x(), y2 = p2.y();

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;
    int err = dx - dy;

    int x = x1, y = y1;

    // Calculate approximate number of points / Vychislyayem priblizitel'noye kolichestvo tochek
    int numPoints = std::max(dx, dy) + 1;
    profile.reserve(numPoints);

    while (true) {
        // Sample intensity at current point / Otshivayem intensivnost' v tekushchey tochke
        if (x >= 0 && x < image.width() && y >= 0 && y < image.height()) {
            profile.append(static_cast<double>(qGray(image.pixel(x, y))));
        } else {
            profile.append(0.0);
        }

        if (x == x2 && y == y2) break;

        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x += sx;
        }
        if (e2 < dx) {
            err += dx;
            y += sy;
        }
    }

    return profile;
}

/*
 * EN: Builds profile along a line between two points with interpolation for subpixel accuracy
 * RU: Stroit profil' vdol' linii mezhdu dvumya tochkami s interpolyatsiyey dlya subpiksel'noy tochnosti
 */
QVector<double> ImageProcessor::buildProfileBetweenPoints(const QPointF& p1, const QPointF& p2,
                                                           const QImage& image, int numSamples)
{
    QVector<double> profile;

    if (image.isNull() || numSamples <= 0) {
        return profile;
    }

    profile.reserve(numSamples);

    double dx = p2.x() - p1.x();
    double dy = p2.y() - p1.y();

    for (int i = 0; i < numSamples; i++) {
        double t = static_cast<double>(i) / (numSamples - 1);
        double x = p1.x() + t * dx;
        double y = p1.y() + t * dy;

        // Bilinear interpolation for subpixel accuracy / Bilineynaya interpolyatsiya dlya subpiksel'noy tochnosti
        int x0 = static_cast<int>(std::floor(x));
        int y0 = static_cast<int>(std::floor(y));
        int x1 = x0 + 1;
        int y1 = y0 + 1;

        double fx = x - x0;
        double fy = y - y0;

        // Bounds checking / Proverka granits
        x0 = std::max(0, std::min(x0, image.width() - 1));
        x1 = std::max(0, std::min(x1, image.width() - 1));
        y0 = std::max(0, std::min(y0, image.height() - 1));
        y1 = std::max(0, std::min(y1, image.height() - 1));

        // Get four surrounding pixel values / Poluchayem znacheniya chetyrekh okruzhayushchikh pikseley
        double v00 = qGray(image.pixel(x0, y0));
        double v01 = qGray(image.pixel(x0, y1));
        double v10 = qGray(image.pixel(x1, y0));
        double v11 = qGray(image.pixel(x1, y1));

        // Bilinear interpolation / Bilineynaya interpolyatsiya
        double v0 = v00 * (1.0 - fx) + v10 * fx;
        double v1 = v01 * (1.0 - fx) + v11 * fx;
        double value = v0 * (1.0 - fy) + v1 * fy;

        profile.append(value);
    }

    return profile;
}

// ============================================================================
// Image Statistics Functions / Funkcii statistiki izobrazheniy
// ============================================================================

/*
 * EN: Converts ImageStatistics to formatted string
 * RU: Preobrazuyet ImageStatistics v formatirovannuyu stroku
 */
QString ImageProcessor::ImageStatistics::toString() const
{
    QString result;
    result += "=== Image Statistics / Statistika izobrazheniya ===\n";
    result += QString("Pixel Count / Kolichestvo pikseley: %1\n").arg(pixelCount);
    result += QString("Sum / Summa: %1\n").arg(sum, 0, 'f', 2);
    result += QString("Mean / Sredneye: %1\n").arg(mean, 0, 'f', 4);
    result += QString("Median / Mediana: %1\n").arg(median, 0, 'f', 4);
    result += QString("Min / Minimum: %1\n").arg(min, 0, 'f', 4);
    result += QString("Max / Maksimum: %1\n").arg(max, 0, 'f', 4);
    result += QString("Range / Razmah: %1\n").arg(max - min, 0, 'f', 4);
    result += QString("Variance / Dispersiya: %1\n").arg(variance, 0, 'f', 4);
    result += QString("Std Deviation / Srednekvadraticheskoye: %1\n").arg(stdDev, 0, 'f', 4);
    result += QString("Skewness / Asimmetriya: %1\n").arg(skewness, 0, 'f', 6);
    result += QString("Kurtosis / Ekstsess: %1\n").arg(kurtosis, 0, 'f', 6);
    result += QString("Entropy / Entropiya: %1\n").arg(entropy, 0, 'f', 6);

    if (edgePixelCount > 0) {
        result += "\n=== Edge Statistics / Statistika granits ===\n";
        result += QString("Edge Pixels / Granichnyye pikseli: %1\n").arg(edgePixelCount);
        result += QString("Edge Ratio / Dolya granits: %1%\n").arg(edgeRatio * 100.0, 0, 'f', 2);
    }

    return result;
}

/*
 * EN: Computes comprehensive statistics for an image
 * RU: Vychislyayet vsestoronyuyu statistiku dlya izobrazheniya
 */
ImageProcessor::ImageStatistics ImageProcessor::computeImageStatistics(
    const QImage& image, const Matrix2D<int>& attributeMatrix)
{
    ImageStatistics stats;

    if (image.isNull()) {
        return stats;
    }

    int width = image.width();
    int height = image.height();
    stats.pixelCount = static_cast<long long>(width) * height;

    // Initialize histogram (256 bins for grayscale) / Initsializiruyem gistogrammu (256 binov dlya seroy shkaly)
    stats.histogram.fill(0, 256);

    // First pass: collect basic stats and histogram / Pervyy prokhod: sbor osnovnoy statistiki i gistogrammy
    double sum = 0.0;
    double sumSq = 0.0;
    stats.min = 255.0;
    stats.max = 0.0;

    // Vectors for median calculation / Vektory dlya vychisleniya mediany
    QVector<quint8> allValues;
    allValues.reserve(stats.pixelCount);

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            int intensity = qGray(image.pixel(i, j));
            allValues.append(static_cast<quint8>(intensity));

            sum += intensity;
            sumSq += static_cast<double>(intensity) * intensity;
            stats.min = std::min(stats.min, static_cast<double>(intensity));
            stats.max = std::max(stats.max, static_cast<double>(intensity));

            // Update histogram / Obnovlyayem gistogrammu
            if (intensity >= 0 && intensity < 256) {
                stats.histogram[intensity]++;
            }
        }
    }

    stats.sum = sum;
    stats.mean = sum / stats.pixelCount;
    stats.variance = (sumSq / stats.pixelCount) - (stats.mean * stats.mean);
    stats.stdDev = std::sqrt(stats.variance);

    // Calculate median / Vychislyayem medianu
    std::sort(allValues.begin(), allValues.end());
    if (allValues.size() % 2 == 0) {
        stats.median = (allValues[allValues.size() / 2 - 1] + allValues[allValues.size() / 2]) / 2.0;
    } else {
        stats.median = allValues[allValues.size() / 2];
    }

    // Calculate skewness and kurtosis (third and fourth moments)
    // Vychislyayem asimmetriyu i ekstsess (tretiy i chetvertyy momenty)
    double m3 = 0.0, m4 = 0.0;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            double diff = qGray(image.pixel(i, j)) - stats.mean;
            m3 += diff * diff * diff;
            m4 += diff * diff * diff * diff;
        }
    }

    if (stats.stdDev > 1e-10) {
        stats.skewness = m3 / (stats.pixelCount * stats.stdDev * stats.stdDev * stats.stdDev);
        stats.kurtosis = m4 / (stats.pixelCount * stats.stdDev * stats.stdDev * stats.stdDev * stats.stdDev) - 3.0;
    }

    // Calculate entropy / Vychislyayem entropiyu
    stats.entropy = 0.0;
    double log2 = std::log(2.0);
    for (int i = 0; i < 256; i++) {
        if (stats.histogram[i] > 0) {
            double p = static_cast<double>(stats.histogram[i]) / stats.pixelCount;
            stats.entropy -= p * std::log(p) / log2;
        }
    }

    // Edge statistics if attribute matrix provided / Statistika granits yesli predostavlena atributnaya matritsa
    if (!attributeMatrix.empty()) {
        for (size_t i = 0; i < attributeMatrix.size(); i++) {
            for (size_t j = 0; j < attributeMatrix[i].size(); j++) {
                if (attributeMatrix[i][j] == static_cast<int>(attribute::isEdge) ||
                    attributeMatrix[i][j] == static_cast<int>(attribute::isSelectedEdge)) {
                    stats.edgePixelCount++;
                }
            }
        }
        stats.edgeRatio = static_cast<double>(stats.edgePixelCount) / stats.pixelCount;
    }

    return stats;
}

/*
 * EN: Computes statistics for a double matrix
 * RU: Vychislyayet statistiku dlya matritsy double
 */
ImageProcessor::ImageStatistics ImageProcessor::computeMatrixStatistics(const Matrix2D<double>& matrix)
{
    ImageStatistics stats;

    if (matrix.empty() || matrix[0].empty()) {
        return stats;
    }

    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    stats.pixelCount = static_cast<long long>(rows) * cols;

    double sum = 0.0;
    double sumSq = 0.0;
    stats.min = std::numeric_limits<double>::max();
    stats.max = -std::numeric_limits<double>::max();

    QVector<double> allValues;
    allValues.reserve(stats.pixelCount);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val = matrix[i][j];
            allValues.append(val);
            sum += val;
            sumSq += val * val;
            stats.min = std::min(stats.min, val);
            stats.max = std::max(stats.max, val);
        }
    }

    stats.sum = sum;
    stats.mean = sum / stats.pixelCount;
    stats.variance = (sumSq / stats.pixelCount) - (stats.mean * stats.mean);
    stats.stdDev = std::sqrt(stats.variance);

    // Calculate median / Vychislyayem medianu
    std::sort(allValues.begin(), allValues.end());
    if (allValues.size() % 2 == 0) {
        stats.median = (allValues[allValues.size() / 2 - 1] + allValues[allValues.size() / 2]) / 2.0;
    } else {
        stats.median = allValues[allValues.size() / 2];
    }

    // Calculate skewness and kurtosis / Vychislyayem asimmetriyu i ekstsess
    double m3 = 0.0, m4 = 0.0;
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double diff = matrix[i][j] - stats.mean;
            m3 += diff * diff * diff;
            m4 += diff * diff * diff * diff;
        }
    }

    if (stats.stdDev > 1e-10) {
        stats.skewness = m3 / (stats.pixelCount * stats.stdDev * stats.stdDev * stats.stdDev);
        stats.kurtosis = m4 / (stats.pixelCount * stats.stdDev * stats.stdDev * stats.stdDev * stats.stdDev) - 3.0;
    }

    return stats;
}

/*
 * EN: Saves statistics to a text file
 * RU: Sohranyayet statistiku v tekstovyy fayl
 */
bool ImageProcessor::saveStatisticsToFile(const ImageStatistics& stats, const QString& filename)
{
    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        return false;
    }

    QTextStream out(&file);
    out << stats.toString();
    file.close();
    return true;
}
















// ============================================================================
// Вспомогательные функции
// ============================================================================

// Билинейная интерполяция значения из матрицы по дробным координатам
double ImageProcessor::interpolateBilinear(const Matrix2D<double>& img, double x, double y)
{
    int rows = static_cast<int>(img.size());
    int cols = static_cast<int>(img[0].size());

    // Clamp to valid range FIRST
    x = std::max(0.0, std::min(x, static_cast<double>(rows - 1)));
    y = std::max(0.0, std::min(y, static_cast<double>(cols - 1)));

    int x0 = static_cast<int>(std::floor(x));
    int y0 = static_cast<int>(std::floor(y));
    int x1 = std::min(x0 + 1, rows - 1);
    int y1 = std::min(y0 + 1, cols - 1);

    x0 = std::max(0, x0);
    y0 = std::max(0, y0);

    double fx = x - x0;
    double fy = y - y0;

    double v00 = img[x0][y0];
    double v10 = img[x1][y0];
    double v01 = img[x0][y1];
    double v11 = img[x1][y1];

    return (1.0 - fx) * (1.0 - fy) * v00 +
           fx  * (1.0 - fy) * v10 +
           (1.0 - fx) *        fy  * v01 +
           fx  *        fy  * v11;
}


// Вычисление остатка между двумя профилями (для градиентного спуска)
double ImageProcessor::computeResidual(const QVector<double>& yy1,
                              const QVector<QPointF>& yP2_full,
                              int XL2, int XR2, int NN)
{
    QVector<QPointF> yyP2 = ImageProcessor::reInterpolateProfile(yP2_full, XL2, XR2, NN);
    QVector<double> yy2;
    yy2.reserve(NN);
    for (int i = 0; i < NN; i++) {
        yy2.push_back(yyP2[i].y());
    }

    double m1 = 0.0, m2 = 0.0;
    for (int i = 0; i < NN; i++) {
        m1 += yy1[i];
        m2 += yy2[i];
    }
    m1 /= NN;
    m2 /= NN;

    double D1 = 0.0, D2 = 0.0;
    for (int i = 0; i < NN; i++) {
        D1 += (yy1[i] - m1) * (yy1[i] - m1);
        D2 += (yy2[i] - m2) * (yy2[i] - m2);
    }
    D1 = std::sqrt(D1 / NN);
    D2 = std::sqrt(D2 / NN);

    if (D1 < 1e-10) D1 = 1.0;
    if (D2 < 1e-10) D2 = 1.0;

    double res = 0.0;
    for (int i = 0; i < NN; i++) {
        double norm1 = (yy1[i] - m1) / D1;
        double norm2 = (yy2[i] - m2) / D2;
        double diff = norm1 - norm2;
        res += diff * diff;
    }
    return res;
}

// ============================================================================
// ОПТИМИЗИРОВАННЫЙ refineSinglePoint
// Теперь принимает ПРЕДВАРИТЕЛЬНО отфильтрованные LoG-изображения
// ============================================================================

RefinementResult ImageProcessor::refineSinglePoint002(int n0, int m0, const RefinementParameters& params)
{
    RefinementResult result;
    result.success = false;
    result.FFF1_0 = 0.0;
    result.FFF1_1 = 0.0;
    result.residual = std::numeric_limits<double>::max();
    result.refinedPosition = QPointF(static_cast<double>(n0), static_cast<double>(m0));

    if (n0 < 0 || n0 >= params.NX || m0 < 0 || m0 >= params.NY) {
        qDebug() << "Error: Point out of bounds";
        return result;
    }

    const int n_sigma = params.n_sigma;
    const int n_myu = params.n_myu;
    const int x0 = n0, y0 = m0;
    const double sigma_myu = std::min(params.sigma1, params.sigma2);
    const double sigma1Sq = params.sigma1 * params.sigma1;

    // Генерация mu
    QVector<double> mu;
    mu.reserve(n_myu);
    double muStep = 2.0 * n_sigma * sigma_myu / n_myu;
    double muStart = -n_sigma * sigma_myu;
    for (int s = 0; s < n_myu; s++) {
        mu.push_back(muStart + s * muStep);
    }

    // ========================================================================
    // ПОСТРОЕНИЕ ПРОФИЛЕЙ (сохранена оригинальная логика)
    // ========================================================================

    QVector<double> prof1, prof2;
    prof1.resize(n_myu);
    prof2.resize(n_myu);

    // Предварительно кэшируем указатели на строки для быстрого доступа
    // (это основная оптимизация для vector<vector<double>>)
    int rows = params.NX;
    int cols = params.NY;

    // Кэшируем exp() для часто используемых значений через LUT
    const int LUT_SIZE = 10000;
    const double MAX_RAD_SQ = 5000.0;
    std::vector<double> expLUT(LUT_SIZE);
    double lutStep = MAX_RAD_SQ / LUT_SIZE;
    double invTwoSigmaSq = 0.5 / sigma1Sq;
    for (int i = 0; i < LUT_SIZE; i++) {
        expLUT[i] = std::exp(-i * lutStep * invTwoSigmaSq);
    }

    // Основной цикл построения профиля
    for (int s = 0; s < n_myu; s++) {
        double prof1Val = 0.0;
        double prof2Val = 0.0;
        double mu_s = mu[s];

        for (int n = 0; n < rows; n++) {
            const std::vector<double>& rowA = params.A[n];
            const std::vector<double>& rowB01 = params.B01[n];
            double dx_base = x0 - n + mu_s * params.ex;

            for (int m = 0; m < cols; m++) {
                double dy = y0 - m + mu_s * params.ey;
                double radSq = dx_base * dx_base + dy * dy;

                // Быстрое вычисление LoG через LUT
                double mult;
                if (radSq < MAX_RAD_SQ) {
                    int lutIdx = static_cast<int>(radSq / lutStep);
                    lutIdx = std::max(0, std::min(lutIdx, LUT_SIZE - 1));
                    mult = (radSq / sigma1Sq - 2.0) * expLUT[lutIdx];
                } else {
                    mult = 0.0;
                }

                prof1Val += rowA[m] * mult;
                prof2Val += rowB01[m] * mult;
            }
        }
        prof1[s] = prof1Val;
        prof2[s] = prof2Val;
    }

    // ========================================================================
    // ВСЁ ОСТАЛЬНОЕ БЕЗ ИЗМЕНЕНИЙ (как в оригинале)
    // ========================================================================

    // Поиск zero-crossings
    App_Stats y2_stats;
    y2_stats.gather_stats(prof2);
    int nmumax = std::max(y2_stats.x_max, y2_stats.x_min);
    int nmumin = std::min(y2_stats.x_max, y2_stats.x_min);

    int nL_Zero = n_myu / 2;
    int nR_Zero = n_myu / 2;
    for (int n = nmumin; n < nmumax && n < n_myu - 1; n++) {
        if (prof1[n] * prof1[n + 1] < 0) nL_Zero = n;
        if (prof2[n] * prof2[n + 1] < 0) nR_Zero = n;
    }

    int N_Zero = (nL_Zero + nR_Zero) / 2;
    const int MAX_COUNTER = 5;
    const int MIN_MARGIN = 20;
    const int SAFE_MARGIN = 35;

    int XL1 = 0, XR1 = 0;
    int decrCounter = 0, incrCounter = 0;

    for (int i = N_Zero + 1; i < n_myu - 1; i++) {
        if (prof2[i] < prof2[i + 1]) incrCounter = std::min(incrCounter + 1, MAX_COUNTER);
        if (prof2[i] > prof2[i + 1]) decrCounter = std::min(decrCounter + 1, MAX_COUNTER);
        if (decrCounter == MAX_COUNTER && incrCounter == MAX_COUNTER) {
            XR1 = std::min(i, n_myu - SAFE_MARGIN);
            break;
        }
    }
    if (XR1 == 0) XR1 = std::min(N_Zero + 50, n_myu - SAFE_MARGIN);

    decrCounter = 0;
    incrCounter = 0;
    for (int i = N_Zero - 1; i > 1; i--) {
        if (prof2[i] < prof2[i + 1]) incrCounter = std::min(incrCounter + 1, MAX_COUNTER);
        if (prof2[i] > prof2[i + 1]) decrCounter = std::min(decrCounter + 1, MAX_COUNTER);
        if (decrCounter == MAX_COUNTER && incrCounter == MAX_COUNTER) {
            XL1 = std::max(i, MIN_MARGIN);
            break;
        }
    }
    if (XL1 == 0) XL1 = std::max(N_Zero - 50, MIN_MARGIN);

    // Подготовка первого профиля
    QVector<QPointF> yP1;
    yP1.reserve(XR1 - XL1);
    for (int i = XL1; i < XR1; i++) {
        yP1.push_back({static_cast<double>(i), prof1[i]});
    }

    const int NN = params.NN;
    QVector<QPointF> yyP1 = reInterpolateProfile(yP1, NN);
    QVector<double> yy1;
    yy1.reserve(NN);
    for (int i = 0; i < NN; i++) {
        yy1.push_back(yyP1[i].y());
    }

    // Поиск оптимальных XL2/XR2 (градиентный спуск вместо полного перебора)
    const int otstup = params.otstup;
    QVector<QPointF> yP2_full;
    yP2_full.reserve(n_myu);
    for (int i = 0; i < n_myu; i++) {
        yP2_full.push_back({static_cast<double>(i), prof2[i]});
    }

    int best_XL2 = XL1, best_XR2 = XR1;
    double minraz = computeResidual(yy1, yP2_full, best_XL2, best_XR2, NN);

    // Градиентный спуск
    for (int step = std::max(1, otstup / 2); step >= 1; step /= 2) {
        bool improved = true;
        while (improved) {
            improved = false;
            for (int dXL = -1; dXL <= 1; dXL++) {
                for (int dXR = -1; dXR <= 1; dXR++) {
                    if (dXL == 0 && dXR == 0) continue;

                    int testXL2 = best_XL2 + dXL * step;
                    int testXR2 = best_XR2 + dXR * step;

                    if (testXL2 < std::max(0, XL1 - otstup) ||
                        testXR2 >= std::min(n_myu, XR1 + otstup) ||
                        testXL2 >= testXR2) continue;

                    double raz = computeResidual(yy1, yP2_full, testXL2, testXR2, NN);

                    if (raz < minraz) {
                        minraz = raz;
                        best_XL2 = testXL2;
                        best_XR2 = testXR2;
                        improved = true;
                    }
                }
            }
        }
    }

    // Финальная обработка
    QVector<QPointF> yyP2_final = reInterpolateProfile(yP2_full, best_XL2, best_XR2, NN);
    QVector<double> yy2_final;
    yy2_final.reserve(NN);
    for (int i = 0; i < NN; i++) {
        yy2_final.push_back(yyP2_final[i].y());
    }

    double m1_f = 0.0, m2_f = 0.0;
    for (int i = 0; i < NN; i++) {
        m1_f += yy1[i];
        m2_f += yy2_final[i];
    }
    m1_f /= NN;
    m2_f /= NN;

    double D1_f = 0.0, D2_f = 0.0;
    for (int i = 0; i < NN; i++) {
        D1_f += (yy1[i] - m1_f) * (yy1[i] - m1_f);
        D2_f += (yy2_final[i] - m2_f) * (yy2_final[i] - m2_f);
    }
    D1_f = std::sqrt(D1_f / NN);
    D2_f = std::sqrt(D2_f / NN);

    if (D1_f < 1e-10) D1_f = 1.0;
    if (D2_f < 1e-10) D2_f = 1.0;

    double MINRAZ_final = 0.0;
    for (int i = 0; i < NN; i++) {
        double norm1 = (yy1[i] - m1_f) / D1_f;
        double norm2 = (yy2_final[i] - m2_f) / D2_f;
        double diff = norm1 - norm2;
        MINRAZ_final += diff * diff;
    }

    double xa = static_cast<double>(XL1);
    double xb = static_cast<double>(XR1);
    double ya = static_cast<double>(best_XL2 - XL1);
    double yb = static_cast<double>(best_XR2 - XR1);

    if (std::abs(yb - ya) < 1e-10) {
        qDebug() << "Error: Division by zero in shift calculation";
        return result;
    }

    double x00 = (xa * yb - xb * ya) / (yb - ya);
    double mu00 = -n_sigma * sigma_myu + x00 * (2.0 * n_sigma * sigma_myu) / n_myu;

    result.FFF1_1 = mu00;
    result.FFF1_0 = (yb - ya) / (xb - xa);
    result.residual = MINRAZ_final;

    double n_new = n0 + params.ex * result.FFF1_1;
    double m_new = m0 + params.ey * result.FFF1_1;
    result.refinedPosition = QPointF(n_new, m_new);
    result.success = true;

    return result;
}



