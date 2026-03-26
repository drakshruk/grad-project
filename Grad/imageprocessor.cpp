#include "imageprocessor.h"

/*
 * EN:
 *  checks if vec contains point p
 * RU:
 *  proveryaet est' li tochka p v vectore vec
 */
bool ImageProcessor::containsVP(QVector<QPoint> &vec, QPoint p)
{
    for(int i = 0; i < vec.size(); i++)
    {
        if(vec[i].x() == p.x() && vec[i].y() == p.y()) return true;
    }
    return false;
}

/*
 * EN:
 *  if pos is the edge pixel
 *  selects all connected pixels into the selectedEdge
 *  and updates iiImAttMat
 * RU:
 *  yesli pos - eto kraevoy piksel' to
 *  vybiraet vse svyazannye piksely v selectedEdge
 *  i obnovlyayet iiImAttMat
 */
void ImageProcessor::selectEdge(QPoint pos, Matrix2D<int> &iiImAttMat, QVector<QPoint> &selectedEdge)
{
    if(!selectedEdge.empty()) selectedEdge.clear();
    if(iiImAttMat[pos.x()][pos.y()] == int(attribute::isEdge))
    {
        selectedEdge.push_back(pos);
        iiImAttMat[pos.x()][pos.y()] = int(attribute::isSelectedEdge);
        size_t prevSize = 0, newSize = 1;
        while(prevSize < newSize)
        {
            prevSize = selectedEdge.size();
            for(int s = 0; s < selectedEdge.size(); s++)
            {
                for(int ii = -1; ii <= 1; ii++)
                {
                    for(int jj = -1; jj <= 1; jj++)
                    {
                        if(iiImAttMat[selectedEdge.at(s).x() + ii][selectedEdge.at(s).y() + jj] == static_cast<int>(attribute::isEdge) &&
                                !containsVP(selectedEdge,QPoint(selectedEdge.at(s).x() + ii, selectedEdge.at(s).y() + jj)))
                        {
                            selectedEdge.push_back(QPoint(selectedEdge.at(s).x() + ii, selectedEdge.at(s).y() + jj));
                            iiImAttMat[selectedEdge.at(s).x() + ii][selectedEdge.at(s).y() + jj] = static_cast<int>(attribute::isSelectedEdge);
                        }
                    }
                }
            }
            newSize = selectedEdge.size();
        }
    }
}

/*
 * EN:
 *  builds profile perpendicular to image edge by taking points along the gradient direction
 * RU:
 *  stroit profil' perpendikulyarno krayu izobrazheniya, berya tochki vdol' napravleniya gradienta
 */
const ProfileResult ImageProcessor::buildProfile001(const ProfileParameters &params, const  Matrix2D<double> &ddImMat)
{
    ProfileResult result;
    double gradientMagnitude = std::sqrt(params.xGradient * params.xGradient + params.yGradient * params.yGradient);
    double normalizedXGrad = params.xGradient / gradientMagnitude;
    double normalizedYGrad = params.yGradient / gradientMagnitude;

    int numSteps = static_cast<int>(params.numSigma * params.sigma);

    for (int step = 0; step < numSteps; ++step) {
        double muS = -params.numSigma * params.sigma +
                         2.0 * params.numSigma * params.sigma * step / numSteps;

        double x = params.x0 + muS * normalizedXGrad;
        double y = params.y0 + muS * normalizedYGrad;

        int xCoord = static_cast<int>(std::round(x));
        int yCoord = static_cast<int>(std::round(y));

        if (xCoord < 0 || xCoord >= (int)ddImMat.size() || yCoord < 0 || yCoord >= (int)ddImMat[0].size()) {
           result.points.push_back(QPointF(static_cast<double>(step), 0.0));
        } else {
           double yProf = ddImMat.at(xCoord).at(yCoord);
           result.points.push_back(QPointF(static_cast<double>(step), yProf));
           result.indices.push_back(QPoint(xCoord, yCoord));
        }
    }
    return result;
}

/*
 * EN:
 *  builds profile perpendicular to image edge by taking points along the gradient direction
 *  and applying laplacian convolution to ddImMat with a kernel of radius convRad

 * RU:
 *  stroit profil' perpendikulyarno krayu izobrazheniya, berya tochki vdol' napravleniya gradienta
 *  i primenyaya svertki laplasa k ddImMat s yadrom radiusa convRad
 */
const ProfileResult ImageProcessor::buildProfile002(const ProfileParameters &params, const Matrix2D<double> &ddImMat, int convRad)
{
    ProfileResult result;

    if (ddImMat.empty() || ddImMat[0].empty()) {
        qDebug() << "Error: Input matrix is empty!";
        return result;
    }

    const int rows = ddImMat.size();
    const int cols = ddImMat[0].size();

    double gradientMagnitude = std::sqrt(params.xGradient * params.xGradient + params.yGradient * params.yGradient);

    if (gradientMagnitude < 1e-10) {
        qDebug() << "Error: Gradient magnitude is zero!";
        return result;
    }

    double normalizedXGrad = params.xGradient / gradientMagnitude;
    double normalizedYGrad = params.yGradient / gradientMagnitude;

//    qDebug() << "1 - Gradient normalized:" << normalizedXGrad << normalizedYGrad;

    int kernelSize = 2 * convRad + 1;
    Matrix2D<double> convCore = getLapl(kernelSize, kernelSize, params.sigma);

    if ((int)convCore.size() != kernelSize || (int)convCore[0].size() != kernelSize) {
//        qDebug() << "Error: Convolution core has wrong size!" << convCore.size() << convCore[0].size();
        return result;
    }

    int numSteps = static_cast<int>(params.numSigma * params.sigma);

    if (numSteps <= 0) {
        numSteps = 50;
//        qDebug() << "Warning: numSteps was" << numSteps << ", using default 50";
    }

    result.points.clear();
    result.indices.clear();
    result.points.reserve(numSteps);
    result.indices.reserve(numSteps);

//    qDebug() << "2 - Starting profile with" << numSteps << "steps";

    for (int step = 0; step < numSteps; ++step) {
        double muS = -params.numSigma * params.sigma +
                     2.0 * params.numSigma * params.sigma * step / numSteps;

        double x = params.x0 + muS * normalizedXGrad;
        double y = params.y0 + muS * normalizedYGrad;

        int xCoord = static_cast<int>(std::round(x));
        int yCoord = static_cast<int>(std::round(y));

//        qDebug() << "Step" << step << ": (" << xCoord << "," << yCoord << ")";

        if (xCoord < 0 || xCoord >= rows || yCoord < 0 || yCoord >= cols) {
            result.points.push_back(QPointF(static_cast<double>(step), 0.0));
            result.indices.push_back(QPoint(xCoord, yCoord));
//            qDebug() << "  Out of bounds -> (0,0)";
        } else {
            double yProf = 0.0;

            for (int i = -convRad; i <= convRad; i++) {
                for (int j = -convRad; j <= convRad; j++) {
                    int kernel_i = i + convRad;
                    int kernel_j = j + convRad;

                    if (kernel_i < 0 || kernel_i >= kernelSize ||
                        kernel_j < 0 || kernel_j >= kernelSize) {
//                        qDebug() << "  Kernel index error:" << kernel_i << kernel_j;
                        continue;
                    }

                    int convX = xCoord + i;
                    int convY = yCoord + j;

                    if(convX >= rows) convX = 2 * (rows - 1) - convX;
                    if(convY >= cols) convY = 2 * (cols - 1) - convY;

                    if (convX >= 0 && convY >= 0) {
                        yProf += convCore[kernel_i][kernel_j] * ddImMat[convX][convY];
                    }
                }
            }

            result.points.push_back(QPointF(static_cast<double>(step), yProf));
            result.indices.push_back(QPoint(xCoord, yCoord));
//            qDebug() << "  Value:" << yProf;
        }
    }

//    qDebug() << "3 - Profile completed with" << result.points.size() << "points";
    return result;
}

/*
 * EN:
 *  cuts the profile from the middle to closest local minimum and maximum
 *  and steps another 'step' pixels to left and right
 * RU:
 *  obrezayet profil' ot seredini do blizhayshih lokal'nogo minimuma i maksimuma
 *  i shagayet yeshche 'step' pikseley vlevo i vpravo
 */
const QVector<QPointF> ImageProcessor::cutProfile(const QVector<QPointF> &profile, int step)
{
    QVector<QPointF> res;
    int left = profile.size()/2, right = profile.size()/2;
    bool leftIncr = true, rightIncr = true, decrLeft = (profile[left-1].y() < profile[left].y());
    while(leftIncr || rightIncr){
        if(left <= 0 || right >= profile.size()-1) break;
        if(decrLeft){
            if(profile[left-1].y() < profile[left].y() && leftIncr){
                left--;
            }
            else if(leftIncr){
                if(left < step) left = 0;
                else left -= step;
                leftIncr = false;
            }
            if(profile[right+1].y() > profile[right].y() && rightIncr){
                right++;
            }
            else if(rightIncr){
                if(right > profile.size() - step - 1) right = profile.size() - 1;
                else right += step;
                rightIncr = false;
            }
        }
        else{
            if(profile[left-1].y() > profile[left].y() && leftIncr){
                left--;
            }
            else if(leftIncr){
                if(left < step) left = 0;
                else left -= step;
                leftIncr = false;
            }
        }
        if(profile[right+1].y() < profile[right].y() && rightIncr){
            right++;
        }
        else if(rightIncr){
            if(right > profile.size() - step - 1) right = profile.size() - 1;
            else right += step;
            rightIncr = false;
        }
    }

    res.clear();
    res.reserve(right - left);

    for(int i = left; i < right; i++){
        res.push_back(profile[i]);
    }
    return res;
}

/*
 * EN:
 *  cuts the profile from left to right
 * RU:
 *  obrezayet profil' sleva napravo
 */
const QVector<QPointF> ImageProcessor::cutProfile(const QVector<QPointF> &profile, int left, int right)
{
    QVector<QPointF> result;
    if(left < 0) left = 0;
    if(right >= profile.size()) right = profile.size() - 1;
    int tmp = left;
    left = left > right ? right : left;
    right = tmp > right ? tmp : right;
    for(int i = left; i <= right; i++){
        if(i < profile.size()){
            result.push_back(profile[i]);
        }
    }
    return result;
}

/*
 * EN:
 *  returns the sum of the square differences of profile1 and profile2 points.y values
 * RU:
 *  vozvrashchaet summu kvadratov raznostey znacheniy y profile1 i profile2
 */
double ImageProcessor::differenceOfTwoProfiles(const QVector<QPointF> &profile1, const QVector<QPointF> &profile2)
{
    if(profile1.size() != profile2.size()) return -1;
    double res = 0.0;
    for(int i = 0; i < profile1.size(); i++){
        res += (profile1[i].y() - profile2[i].y()) * (profile1[i].y() - profile2[i].y());
    }
    return res;
}


/*
 * EN:
 *  reinterpolates profile from left bound to right bound
 *  to the new length
 *  using linear interpolation
 * RU:
 *  pereinterpolirayet profil' do novoy dliny
 *  ispol'zuya lineynuyu interpolyatsiyu
 */
const QVector<QPointF> ImageProcessor::reInterpolateProfile(const QVector<QPointF> &profile, int newLen)
{
    return reInterpolateProfile(profile, 0.0, profile.size() - 1, newLen);
};

/*
 * EN:
 *  reinterpolates profile from left bound to right bound
 *  to the new length
 *  using linear interpolation
 * RU:
 *  pereinterpolirayet profil' nachinaya s levoi granitsy do pravoi do novoy dliny
 *  ispol'zuya lineynuyu interpolyatsiyu
 */
const QVector<QPointF> ImageProcessor::reInterpolateProfile(const QVector<QPointF> &profile, double leftBound, double rightBound, int newLen)
{
    QVector<QPointF> newProfile;

    if (profile.size() < 2 || newLen <= 0) {
        qDebug() << "Error: Invalid input in reInterpolateProfile";
        return profile;
    }

    leftBound = std::max(leftBound, 0.0);
    rightBound = std::min(rightBound, static_cast<double>(profile.size() - 1));

    if (leftBound >= rightBound) {
        qDebug() << "Error: Invalid bounds in reInterpolateProfile";
        return profile;
    }

    newProfile.reserve(newLen);

    double step = (rightBound - leftBound) / (newLen - 1);

    for (int i = 0; i < newLen; i++) {
        double currentPos = leftBound + i * step;

        int idx1 = static_cast<int>(std::floor(currentPos));
        int idx2 = static_cast<int>(std::ceil(currentPos));

        if (idx1 < 0) idx1 = 0;
        if (idx2 >= profile.size()) idx2 = profile.size() - 1;
        if (idx1 == idx2) {
            newProfile.append(QPointF(i, profile[idx1].y()));
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

/*
 * EN:
 *  converts image to gray color
 *  and puts it into the result matrix
 * RU:
 *  preobrazuyet izobrazhenie v seriy tsvet
 *  i pomeshchaet v rezul'tiruyushchuyu matricu
 */
const Matrix2D<double> ImageProcessor::fromGrayImage(const QImage &image)
{
    Matrix2D<double> result(image.width(), std::vector<double>(image.height()));
    for(int i = 0; i < image.width(); i++)
    {
        for(int j = 0; j < image.height(); j++)
        {
            result[i][j] = 1.0 * qGray(image.pixel(i,j));
        }
    }
    return result;
}

/*
 * EN:
 *  constructs image from ddImageMat
 *  with gray pixels
 * RU:
 *  sozdayet izobrazhenie iz ddImageMat
 *  s serimi pikselyami
 */
const QImage ImageProcessor::toGrayImage(const Matrix2D<double> &ddImageMat)
{
    QImage resIm(ddImageMat.size(), ddImageMat[0].size(), QImage::Format_ARGB32);

    double dMin = min(ddImageMat);
    double dMax = max(ddImageMat);

    Matrix2D<double> ddTmpMat = ddImageMat;

    if(dMin < 0)
    {
        for(size_t i = 0; i < ddImageMat.size(); i++)
        {
            for(size_t j = 0; j < ddImageMat[0].size(); j++)
            {
                ddTmpMat[i][j] += dMin;
                ddTmpMat[i][j] *= (255 / (dMax + dMin));
                resIm.setPixel(i,j,qRgb(static_cast<int>(ddTmpMat[i][j]), static_cast<int>(ddTmpMat[i][j]),static_cast<int>(ddTmpMat[i][j])));
            }
        }
    }
    if(dMin >= 0)
    {
        for(size_t i = 0; i < ddImageMat.size(); i++)
        {
            for(size_t j = 0; j < ddImageMat[0].size(); j++)
            {
                ddTmpMat[i][j] *= (255 / dMax);
                resIm.setPixel(i,j,qRgb(static_cast<int>(ddTmpMat[i][j]), static_cast<int>(ddTmpMat[i][j]),static_cast<int>(ddTmpMat[i][j])));
            }
        }
    }
    return resIm;
}

/*
 * EN:
 *  constructs image from ddImageMat
 *  turns positive pixels to red
 *  and negative pixels to blue
 * RU:
 *  sozdayet izobrazhenie iz ddImageMat
 *  preobrazuyet polozhitel'nye piksely v krasnyy
 *  a otritsatel'nye piksely v siniy
 */
const QImage ImageProcessor::toBlueRedImage(const Matrix2D<double> &ddImageMat, double dRedMax, double dBlueMax)
{

    QImage resIm(ddImageMat.size(), ddImageMat[0].size(), QImage::Format_ARGB32);

    double dMin = min(ddImageMat);
    double dMax = max(ddImageMat);
    Matrix2D<double> ddTmpMat = ddImageMat;

    for(size_t i = 0; i < ddImageMat.size(); i++)
    {
        for(size_t j = 0; j < ddImageMat.size(); j++)
        {
            if(ddImageMat[i][j] < 0)
            {
                ddTmpMat[i][j] *= (dBlueMax/dMin);
                resIm.setPixel(i,j,qRgb(0, 0,static_cast<int>(ddTmpMat[i][j])));
            }
            else if(ddImageMat[i][j] > 0)
            {
                ddTmpMat[i][j] *= (dRedMax/dMax);
                resIm.setPixel(i,j,qRgb(static_cast<int>(ddTmpMat[i][j]), 0, 0));
            }
            else resIm.setPixel(i,j,qRgb(0,0,0));
        }
    }
    return resIm;
}


/*
 * EN:
 *  updates Image with (255,255,255) pixels if iiImageAttMat = isEdge
 *  and (0,255,255) pixels if iiImageAttMat = isSelectedEdge
 * RU:
 *  obnovlyayet Image s pikselyami (255,255,255) yesli iiImageAttMat = isEdge
 *  i s pikselyami (0,255,255) yesli iiImageAttMat = isSelectedEdge
 */
const QImage ImageProcessor::combineImageWithEdge(const QImage &Image, const Matrix2D<int> &iiImageAttMat)
{
    QImage resIm(Image.width(), Image.height(), QImage::Format_ARGB32);
    for(int i = 0; i < Image.width(); i++){
        for(int j = 0; j < Image.height(); j++){
            if(iiImageAttMat[i][j] == static_cast<int>(attribute::isEdge)) resIm.setPixel(i,j,qRgb(255,255,255));
            else if(iiImageAttMat[i][j] == static_cast<int>(attribute::isSelectedEdge)) resIm.setPixel(i,j,qRgb(0,255,255));
            else resIm.setPixel(i,j, Image.pixel(i,j));
        }
    }
    return resIm;
}

/*
 * EN:
 *  constructs QImage with (255,255,255) pixels if iiImageAttMat = isEdge
 *  and (0,255,255) pixels if iiImageAttMat = isSelectedEdge
 *  and ddImageMat values if not edge
 * RU:
 *  sozdayet QImage s pikselyami (255,255,255) yesli iiImageAttMat = isEdge
 *  i s pikselyami (0,255,255) yesli iiImageAttMat = isSelectedEdge
 *  i znacheniyami ddImageMat yesli ne kray
 */
const QImage ImageProcessor::combineImageWithEdge(const Matrix2D<double> &ddImageMat, const Matrix2D<int> &iiImageAttMat)
{
    QImage resIm(ddImageMat.size(), ddImageMat[0].size(), QImage::Format_ARGB32);
    for(size_t i = 0; i < ddImageMat.size(); i++){
        for(size_t j = 0; j < ddImageMat[0].size(); j++){
            if(iiImageAttMat[i][j] == static_cast<int>(attribute::isEdge)) resIm.setPixel(i,j,qRgb(255,255,255));
            else if(iiImageAttMat[i][j] == static_cast<int>(attribute::isSelectedEdge)) resIm.setPixel(i,j,qRgb(0,255,255));
            else resIm.setPixel(i,j, qRgb(ddImageMat[i][j], ddImageMat[i][j], ddImageMat[i][j]));
        }
    }
    return resIm;
}

/*
 * EN:
 *  updates image im1 with profileIndices to a color (0, 255, 255)
 * RU:
 *  obnovlyayet izobrazhenie im1 s profileIndices tsvetom (0, 255, 255)
 */
const QImage ImageProcessor::combineImageWithProfile(const QImage &im1, const QVector<QPoint> &profileIndices)
{
    QImage resIm(im1.width(), im1.height(), QImage::Format_ARGB32);
    for(int i = 0; i < profileIndices.size(); i++){
        resIm.setPixel(profileIndices[i].x(),profileIndices[i].y(),qRgb(0,255,255));
    }
    return resIm;
}

/*
 * EN:
 *  convolutes image with ddConvCore kernel
 * RU:
 *  svertvyayet izobrazhenie s yadrom ddConvCore
 */
const QImage ImageProcessor::convImage(const QImage &image, const Matrix2D<double> &ddConvCore)
{
    QImage resIm = image;
    int iRad = ddConvCore.size();
    double lSum = 0.0;
    for(int ii = 0; ii < iRad; ii++) {
        for(int jj = 0; jj < iRad; jj++) {
            lSum += ddConvCore[ii][jj];
        }
    }
    if (lSum <= 0) lSum = 1.0;

    for(int i = 0; i < image.width(); i++)
    {
        for(int j = 0; j < image.height(); j++)
        {
            double rSum = 0, gSum = 0, bSum = 0;
            for(int ii = 0; ii < iRad; ii++)
            {
                int it = abs(i - iRad/2 + ii);
                if( it >= image.width())
                {
                    it = 2*(image.width() - 1) - it;
                }
                for(int jj = 0; jj < iRad; jj++)
                {
                    int jt = abs(j - iRad/2 + jj);

                    if( jt >= image.height())
                    {
                        jt = 2*(image.height() - 1) - jt;
                    }
                    rSum += static_cast<double>(qRed(image.pixel(it,jt))) * ddConvCore[ii][jj];
                    gSum += static_cast<double>(qGreen(image.pixel(it,jt))) * ddConvCore[ii][jj];
                    bSum += static_cast<double>(qBlue(image.pixel(it,jt))) * ddConvCore[ii][jj];
                }
            }
            rSum /= lSum;
            if (rSum < 0) rSum = 0;
            if (rSum > 255) rSum = 255;
            gSum /= lSum;
            if (gSum < 0) gSum = 0;
            if (gSum > 255) gSum = 255;
            bSum /= lSum;
            if (bSum < 0) bSum = 0;
            if (bSum > 255) bSum = 255;
            resIm.setPixel(i,j,qRgb(rSum,gSum,bSum));
        }
    }
    return resIm;
}

/*
 * EN:
 *  convolutes ddImageMat matrix with ddConvCore kernel
 * RU:
 *  svertvyayet matricu ddImageMat s yadrom ddConvCore
 */
const Matrix2D<double> ImageProcessor::convMat(const Matrix2D<double> &ddImageMat, const Matrix2D<double> &ddConvCore)
{
    int rows = ddImageMat.size();
    int cols = ddImageMat[0].size();
    int kernelSize = ddConvCore.size();
    int kernelRadius = kernelSize / 2;

    Matrix2D<double> ddResMat(rows, std::vector<double>(cols, 0.0));

    double kernelSum = 0.0;
    for(int i = 0; i < kernelSize; i++) {
        for(int j = 0; j < kernelSize; j++) {
            kernelSum += ddConvCore[i][j];
        }
    }
    if (kernelSum <= 0) kernelSum = 1.0;

    std::vector<int> xIndices(rows + 2*kernelRadius);
    std::vector<int> yIndices(cols + 2*kernelRadius);

    for(int i = -kernelRadius; i < rows + kernelRadius; i++) {
        int idx = i;
        if(i < 0) idx = -i;
        else if(i >= rows) idx = 2*(rows - 1) - i;
        xIndices[i + kernelRadius] = idx;
    }

    for(int j = -kernelRadius; j < cols + kernelRadius; j++) {
        int idx = j;
        if(j < 0) idx = -j;
        else if(j >= cols) idx = 2*(cols - 1) - j;
        yIndices[j + kernelRadius] = idx;
    }

    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            double matSum = 0.0;

            for(int ki = 0; ki < kernelSize; ki++) {
                int imgX = xIndices[i - kernelRadius + ki + kernelRadius];

                for(int kj = 0; kj < kernelSize; kj++) {
                    int imgY = yIndices[j - kernelRadius + kj + kernelRadius];

                    matSum += ddImageMat[imgX][imgY] * ddConvCore[ki][kj];
                }
            }

            ddResMat[i][j] = matSum / kernelSum;
        }
    }

    return ddResMat;
}

/*
 * EN:
 *  fills the matrix with gaussian distibution function values
 * RU:
 *  zapolnyayet matricu funktsiey raspredeleniya Gaussa
 */
const Matrix2D<double> ImageProcessor::getGauss(int iXsize, int iYsize, double dSigma)
{
    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize, 0.0));
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    double sum = 0.;
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = abs(xCenter - i);
            int dy = abs(yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = 1.0 / (dSigma * dSigma * 2.0 * pi) * exp((-1.0 / 2.0) * (Rad * Rad) / (dSigma * dSigma));
            sum += ddRes[i][j];
        }
    }
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            ddRes[i][j] /= sum;
        }
    }
    qDebug() << "gauss kernel sum = " << sum << "\n";
    return ddRes;
}

/*
 * EN:
 *  fills the matrix with gaussian distribution function values
 *  with coefficient dx = x distance from the center
 * RU:
 *  zapolnyayet matricu znacheniyami funktsii raspredeleniya Gaussa
 *  s koeffitsientom dx = rasstoyanie x ot tsentra
 */
const Matrix2D<double> ImageProcessor::getXGradCore(int iXsize, int iYsize, double dSigma)
{
    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize, 0.0));
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    double sum = 0.0;
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = (i - xCenter);
            int dy = (j - yCenter);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = dx / dSigma / dSigma * exp((-1.0 / 2.0) * (Rad * Rad) / (dSigma * dSigma));
            sum += ddRes[i][j];
        }
    }
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
//            ddRes[i][j] /= sum;
        }
    }
    qDebug() << "x grad kernel sum = " << sum << "\n";
    return ddRes;
}

/*
 * EN:
 *  fills the matrix with gaussian distribution function values
 *  with coefficient dy = y distance from the center
* RU:
 *  zapolnyayet matricu znacheniyami funktsii raspredeleniya Gaussa
 *  s koeffitsientom dy = rasstoyanie y ot tsentra
 */
const Matrix2D<double> ImageProcessor::getYGradCore(int iXsize, int iYsize, double dSigma)
{
    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize, 0.0));
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    double sum = 0.0;
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = (i - xCenter);
            int dy = (j - yCenter);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = dy / dSigma / dSigma * exp((-1.0 / 2.0) * (Rad * Rad) / (dSigma * dSigma));
            sum += ddRes[i][j];
        }
    }
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
//            ddRes[i][j] /= sum;
        }
    }
    qDebug() << "y grad kernel sum = " << sum << "\n";
    return ddRes;
}

/*
 * EN:
 *  fills the matrix with laplacian distributions function values
 * RU:
 *  zapolnyayet matricu znacheniyami funktsii raspredeleniya Laplasa
 */
const Matrix2D<double> ImageProcessor::getLapl(int iXsize, int iYsize, double dSigma)
{
    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize, 0.0));
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    double sum = 0.0;
//    double normCoef = 2. / sqrt(3. * dSigma) / sqrt(sqrt(pi));
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = abs(xCenter - i);
            int dy = abs(yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
//            ddRes[i][j] = normCoef * (1. - Rad * Rad / dSigma / dSigma) * exp( -1. * Rad * Rad / 2 / dSigma / dSigma);
            ddRes[i][j] = (Rad * Rad / dSigma / dSigma - 2.0) * exp( -1. * Rad * Rad / 2 / dSigma / dSigma);
            sum += ddRes[i][j];
        }
    }
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
//            ddRes[i][j] /= sum;
        }
    }
    qDebug() << "lapl kernel sum = " << sum << "\n";
    return ddRes;
}

/*
 * EN:
 *  finds all pixels with both positive and negative neighbours in iRad radius
 *  also updates the attribute matrix
 * RU:
 *  nahodit vse piksely s polozhitel'nymi i otritsatel'nymi sosedyami v radiuse iRad
 *  takzhe obnovlyayet matricu atributov
 */
const Matrix2D<double> ImageProcessor::findEdges(const Matrix2D<double> &ddImageMat, int iRad, Matrix2D<int> &iiImageAttMat)
{
    int rows = ddImageMat.size();
    int cols = ddImageMat[0].size();
    Matrix2D<double> ddImageEdgeMat(rows, std::vector<double>(cols, 0.0));

    bool bPos, bNeg;
//    qDebug() << 1 << rows << " " << cols;
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            bPos = false; bNeg = false;

//            qDebug() << i << " " << j;

//            if (i >= 56 && j >= 121) {
//                qDebug() << "=== DEBUG Pixel (" << i << "," << j << ") ===";
//            }

            for(int ii = 0; ii < iRad; ii++)
            {
                for(int jj = 0; jj < iRad; jj++)
                {
                    int it = abs(i - iRad/2 + ii);
                    int jt = abs(j - iRad/2 + jj);
                    if( it >= rows)
                    {
                        it = 2*(rows - 1) - it;
                    }
                    if( jt >= cols)
                    {
                        jt = 2*(cols - 1) - jt;
                    }
                    if (it < 0 || it >= rows || jt < 0 || jt >= cols) {
                         qDebug() << "CRITICAL: Mirroring failed! it:" << it << "jt:" << jt
                                  << "for pixel (" << i << "," << j << ") with offset (" << ii - iRad/2 << "," << jj - iRad/2 << ")";
                         continue;
                     }

//                    qDebug() << it << " " << jt;
//                    qDebug() << ddImageMat[it][jt];

//                    if (i == 56 && j == 121) {
//                          qDebug() << "  Offset (" << ii - iRad/2 << "," << jj - iRad/2 << ") -> (" << it << "," << jt
//                                   << ") value:" << ddImageMat[it][jt];
//                      }

                    if(ddImageMat[it][jt] > 0.) bPos = true;
                    else if(ddImageMat[it][jt] < 0.) bNeg = true;
                }
            }


            if(bPos && bNeg)
            {
//                qDebug() << "=== Trying to update atribute ===";
                iiImageAttMat[i][j] = static_cast<int>(attribute::isEdge);
//                qDebug() << "=== Attribute at: (" << i << "," << j <<"): " << iiImageAttMat[i][j] << " ===";
                ddImageEdgeMat[i][j] = 255.;
            }
            else {
                iiImageAttMat[i][j] = 0.;
                ddImageEdgeMat[i][j] = 0.;
            }
//            if (i == 56 && j == 121) {
//                 qDebug() << "Result: bPos =" << bPos << "bNeg =" << bNeg << "-> value:" << ddImageEdgeMat[i][j];
//                 qDebug() << "=== END DEBUG ===";
//             }
        }
    }
    return ddImageEdgeMat;
}

/*
 * EN:
 *  finds all pixels with both positive and negative neighbours in iRad radius
 * RU:
 *  nahodit vse piksely s polozhitel'nymi i otritsatel'nymi sosedyami v raduse iRad
 */
const Matrix2D<double> ImageProcessor::findEdges(const Matrix2D<double> &ddImageMat, int iRad)
{
    int rows = ddImageMat.size();
    int cols = ddImageMat[0].size();
    Matrix2D<double> ddImageEdgeMat(rows, std::vector<double>(cols, 0.0));
    bool bPos, bNeg;

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            bPos = false; bNeg = false;
            for(int ii = 0; ii < iRad; ii++)
            {
                for(int jj = 0; jj < iRad; jj++)
                {
                    int it = i + (ii - iRad/2);
                    int jt = j + (jj - iRad/2);

                    if (it < 0) it = -it;
                    else if (it >= rows) it = 2 * (rows - 1) - it;

                    if (jt < 0) jt = -jt;
                    else if (jt >= cols) jt = 2 * (cols - 1) - jt;

                    if (it < 0 || it >= rows || jt < 0 || jt >= cols) {
                         qDebug() << "CRITICAL: Mirroring failed! it:" << it << "jt:" << jt
                                  << "for pixel (" << i << "," << j << ") with offset (" << ii - iRad/2 << "," << jj - iRad/2 << ")";
                         continue;
                     }

                    if(ddImageMat[it][jt] > 0) bPos = true;
                    else if(ddImageMat[it][jt] < 0) bNeg = true;
                }
            }
            if(bPos && bNeg)
            {
                ddImageEdgeMat[i][j] = 255.;
            }
            else ddImageEdgeMat[i][j] = 0.;
        }
    }
    return ddImageEdgeMat;
}

/*
 * EN:
 *  finds all pixels with both positive and negative neighbours in iRad radius
 *  if the said pixel > 0 then it is set to red
 *  else it is set to blue
 *  also updates the attribute matrix
 * RU:
 *  nahodit vse piksely s polozhitel'nymi i otritsatel'nymi sosedyami v raduse iRad
 *  yesli takoy piksel' > 0 to on ustanavlivayetsya v krasnyy
 *  inache on ustanavlivayetsya v siniy
 *  takzhe obnovlyayet matricu atributov
 */
const Matrix2D<double> ImageProcessor::findEdgesRB(const Matrix2D<double> &ddImageMat, int iRad, Matrix2D<int> &iiImageAttMat)
{
    int rows = ddImageMat.size();
    int cols = ddImageMat[0].size();
    Matrix2D<double> ddImageEdgeMat(rows, std::vector<double>(cols, 0.0));
    bool bPos, bNeg;
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            bPos = false; bNeg = false;
            for(int ii = 0; ii < iRad; ii++)
            {
                for(int jj = 0; jj < iRad; jj++)
                {
                    int it = i + (ii - iRad/2);
                    int jt = j + (jj - iRad/2);

                    if (it < 0) it = -it;
                    else if (it >= rows) it = 2 * (rows - 1) - it;

                    if (jt < 0) jt = -jt;
                    else if (jt >= cols) jt = 2 * (cols - 1) - jt;

                    if (it < 0 || it >= rows || jt < 0 || jt >= cols) {
                         qDebug() << "CRITICAL: Mirroring failed! it:" << it << "jt:" << jt
                                  << "for pixel (" << i << "," << j << ") with offset (" << ii - iRad/2 << "," << jj - iRad/2 << ")";
                         continue;
                     }

                    if(ddImageMat[it][jt] > 0) bPos = true;
                    else if(ddImageMat[it][jt] < 0) bNeg = true;
                }
            }
            if(bPos && bNeg)
            {
                if(ddImageMat[i][j] > 0) ddImageEdgeMat[i][j] = 1;
                if(ddImageMat[i][j] <= 0) ddImageEdgeMat[i][j] = -1;
                iiImageAttMat[i][j] = static_cast<int>(attribute::isEdge);
            }
            else {
                ddImageEdgeMat[i][j] = 0;
                iiImageAttMat[i][j] = 0;
            }
        }
    }
    return ddImageEdgeMat;
}

/*
 * EN:
 *  gaussian edge detection but also upgrades image attribute matrix
 * RU:
 *  opredelenie krayev po Gaussu, takzhe obnovlyayet matricu atributov izobrazheniya
 */
const QImage ImageProcessor::gaussianEdgeDetection(const Matrix2D<double> &ddImageMat, double dSigma, int iRad, Matrix2D<int> &iiImageAttMat, Matrix2D<double> &ddMatForProfile)
{
    Q_UNUSED(ddMatForProfile);
    Matrix2D<double> mat1;
    Matrix2D<double> mat2;

    mat1 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma));
    mat2 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma*1.6));
    mat1 = elementWiseOperation(mat1, mat2, MatrixLambdas::Subtract<double>{});

    mat1 = findEdges(mat1, 3, iiImageAttMat);
    QImage image = toGrayImage(mat1);
    return image;
}

/*
 * EN:
 *  edge detection using gradient method
 *  convolute image with two gaussians with two different sigma coefficients and substract them
 *  then finds edges using findEdges method
 *  updates iiImageAttMat and saves ddMatForProfile for later use
 * RU:
 *  opredelenie krayev s ispol'zovaniem gradientnogo metoda
 *  svertka izobrazheniya s dvumya funktsiyami Gaussa s razlichnymi sigma-koeffitsientami i ikh vychitanie
 *  zatem nahodit kraya s ispol'zovaniem metoda findEdges
 *  obnovlyayet iiImageAttMat i sohranyaet ddMatForProfile dlya dal'neyshego ispol'zovaniya
 */
const QImage ImageProcessor::gaussianEdgeDetection(const Matrix2D<double> &ddImageMat, double dSigma, int iRad)
{
    Matrix2D<double> mat1;
    Matrix2D<double> mat2;

    mat2 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma));
    mat1 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma*2.0));
    mat1 = elementWiseOperation(mat1, mat2, MatrixLambdas::Subtract<double>{});
    mat1 = findEdges(mat1, 3);
    QImage image = toGrayImage(mat1);
    return image;
}

/*
 * EN:
 *
 * RU:
 *
 */
const QImage ImageProcessor::sampleTwoHollowsBig()
{
    int tmp = 0;
    QImage res(400,550,QImage::Format_ARGB32);
    for(int i = 0; i < 400; i++){
        for(int j = 0; j < 550; j++){
            res.setPixel(i,j,qRgb(0,0,0));
            if(sqrt((i-200.)*(i-200.) + (j-200.)*(j-200.)) < 100.){
                tmp = sqrt(100.*100.-(i-200.)*(i-200.)-(j-200.)*(j-200.));
                res.setPixel(i,j,qRgb(tmp,tmp,tmp));
            }
            if(sqrt((i-200.)*(i-200.) + (j-375.)*(j-375.)) < 75.){
               tmp = sqrt(75.*75.-(i-200.)*(i-200.)-(j-375.)*(j-375.));
               res.setPixel(i,j,qRgb(tmp,tmp,tmp));
            }
        }
    }
    return res;
}

/*
 * EN:
 *
 * RU:
 *
 */
const QImage ImageProcessor::sampleTwoHollows()
{
    int tmp = 0;
    QImage res(200,200,QImage::Format_ARGB32);
    for(int i = 0; i < 200; i++){
        for(int j = 0; j < 200; j++){
            res.setPixel(i,j,qRgb(0,0,0));
            if(sqrt((i-100.)*(i-100.) + (j-135.)*(j-135.)) < 40.){
                tmp = sqrt(40.*40.-(i-100.)*(i-100.)-(j-135.)*(j-135.));
                res.setPixel(i,j,qRgb(tmp,tmp,tmp));
            }
            if(sqrt((i-100.)*(i-100.) + (j-65.)*(j-65.)) < 30.){
               tmp = sqrt(30.*30.-(i-100.)*(i-100.)-(j-65.)*(j-65.));
               res.setPixel(i,j,qRgb(tmp,tmp,tmp));
            }
        }
    }
    return res;
}

/*
 * EN:
 *
 * RU:
 *
 */
double ImageProcessor::calculateResidualOfProfiles(const QVector<QPointF> &profile1, const QVector<QPointF> &profile2)
{
    if (profile1.size() != profile2.size() || profile1.isEmpty()) {
        return std::numeric_limits<double>::max();
    }


    long double sumSquared = 0.0;
    for (int i = 0; i < profile1.size(); ++i) {
        long double diff = profile1[i].y() - profile2[i].y();
        sumSquared += diff * diff;
    }

    return sqrt(sumSquared);

}

/*
 * EN:
 *
 * RU:
 *
 */
RefinementResult ImageProcessor::refineSinglePoint(int n0, int m0, RefinementParameters params)
{
    RefinementResult result;
    /*
     *   QPointF refinedPosition;    // Refined position / Utochnennaya poziciya
     *   double FFF1_1;              // Shift value / Velichina sdviga
     *   double FFF1_0;              // Scale value / Velichina masshtaba
     *   double residual;            // Final residual / Finalnaya nevyazka
     *   bool success;
    */

    // Initialize result with default values / Inicializaciya rezultata znacheniyami po-umolchaniyu
    result.success = false;
    result.FFF1_0 = 0.0;
    result.FFF1_1 = 0.0;
    result.residual = std::numeric_limits<double>::max();
    result.refinedPosition = QPointF(n0, m0);

    // Check if point is within image bounds / Proverka, nahoditsya li tochka v predelah izobrazheniya
    if (n0 < 0 || n0 >= params.NX || m0 < 0 || m0 >= params.NY) {
        qDebug() << "Error: Point out of bounds / Oshibka: Tochka vne granic";
        return result;
    }

    QVector<double> prof1, prof2;
    int n_sigma = params.n_sigma;           // Number of sigma steps / Kolichestvo shagov sigma
    int n_myu = params.n_myu;                // Number of mu points / Kolichestvo tochek mu
    int x0 = n0, y0 = m0;                    // Starting point / Nachalnaya tochka
    double sigma_myu = min(params.sigma1, params.sigma2);  // Sigma for profile / Sigma dlya profilya
    QVector<double> mu;

    // Reserve memory for vectors / Rezerviruem pamyat dlya vektorov
    mu.reserve(n_myu);
    prof1.reserve(n_myu);
    prof2.reserve(n_myu);

    // Generate mu values / Generiruem znacheniya mu
    for(int s = 0; s < n_myu; s++) {
        double val = -n_sigma*sigma_myu + s*(2*n_sigma*sigma_myu)/n_myu;
        mu.push_back(val);
    }

    // Build profiles along gradient direction / Stroim profili vdol napravleniya gradienta
    for(int s = 0; s < n_myu; s++) {
        double prof1Val = 0.0;
        double prof2Val = 0.0;

        for(int n = 0; n < params.NX; n++) {
            for(int m = 0; m < params.NY; m++) {
                // Calculate distance from point along gradient / Vychislyaem rasstoyanie ot tochki vdol gradienta
                double dx = (x0 - n + mu[s]*params.ex);
                double dy = (y0 - m + mu[s]*params.ey);
                double rad = dx*dx + dy*dy;

                // Laplacian of Gaussian kernel / Yadro Laplaciana Gaussa
                double mult = (rad/params.sigma1/params.sigma1 - 2.0) * exp(-0.5*rad/params.sigma1/params.sigma1);

                prof1Val += params.A[n][m] * mult;
                prof2Val += params.B01[n][m] * mult;
            }
        }
        prof1.push_back(prof1Val);
        prof2.push_back(prof2Val);
    }

    QVector<double> y1, y2;  // Only need y values for profiles / Nuzhny tolko znacheniya y dlya profiley
    for(int i = 0; i < n_myu; i++) {
        y1.push_back(prof1[i]);
        y2.push_back(prof2[i]);
    }

    // Find zero crossings / Nahodim perehody cherez nol
    App_Stats y2_stats;
    y2_stats.gather_stats(y2);
    int nmumax = max(y2_stats.x_max, y2_stats.x_min);
    int nmumin = min(y2_stats.x_max, y2_stats.x_min);
    int N0 = n_myu;

    // Add bounds checking for array access / Dobavlyaem proverku granic pri dostupe k massivu
    int nL_Zero = N0/2, nR_Zero = N0/2;
    for(int n = nmumin; n < nmumax && n < n_myu - 1; n++) {
        if(prof1[n]*prof1[n+1] < 0) nL_Zero = n;
        if(prof2[n]*prof2[n+1] < 0) nR_Zero = n;
    }

    int N_Zero = (nL_Zero + nR_Zero)/2;

    // Find profile bounds by counting increases/decreases / Nahodim granicy profilya, schitaya uvelicheniya/umensheniya
    // Using second profile (more blurred) for better stability / Ispolzuem vtoroy profil (bolee razmytyy) dlya luchshey stabilnosti
    int decrCounter = 0, incrCounter = 0;
    int XL1 = 0, XR1 = 0;

    // Search right bound / Poisk pravoy granicy
    for(int i = N_Zero + 1; i < n_myu - 1; i++) {
        if(prof2[i] < prof2[i+1]) {
            incrCounter++;
            incrCounter = min(incrCounter, 5);  // Limit to 5 / Ogranichivaem do 5
        }
        if(prof2[i] > prof2[i+1]) {
            decrCounter++;
            decrCounter = min(decrCounter, 5);
        }
        if(decrCounter == 5 && incrCounter == 5) {
            XR1 = min(i, n_myu - 35);
            break;
        }
    }

    // Add default value if loop didn't find bounds / Dobavlyaem znachenie po-umolchaniyu esli cikl ne nashel granicy
    if (XR1 == 0) XR1 = min(N_Zero + 50, n_myu - 35);

    // Search left bound / Poisk levoy granicy
    incrCounter = 0; decrCounter = 0;
    for(int i = N_Zero - 1; i > 1; i--) {
        if(prof2[i] < prof2[i+1]) {
            incrCounter++;
            incrCounter = min(incrCounter, 5);
        }
        if(prof2[i] > prof2[i+1]) {
            decrCounter++;
            decrCounter = min(decrCounter, 5);
        }
        if(decrCounter == 5 && incrCounter == 5) {
            XL1 = max(i, 20);
            break;
        }
    }

    // Add default value if loop didn't find bounds / Dobavlyaem znachenie po-umolchaniyu esli cikl ne nashel granicy
    if (XL1 == 0) XL1 = max(N_Zero - 50, 20);

    // Prepare first profile for interpolation / Podgotavlivaem pervyy profil dlya interpolyacii
    QVector<QPointF> yP1;
    yP1.reserve(XR1 - XL1);
    for(int i = XL1; i < XR1; i++) {
        yP1.push_back({ 1. * i, prof1[i]} );
    }

    int NN = params.NN;  // New profile length / Novaya dlina profilya

    // Reinterpolate first profile / Pereinterpoliruem pervyy profil
    QVector<QPointF> yyP1 = ImageProcessor::reInterpolateProfile(yP1, NN);
    QVector<double> yy1;
    yy1.reserve(NN);
    for(int s1 = 0; s1 < NN; s1++) {
        yy1.push_back(yyP1[s1].y());
    }

    // Search for optimal second profile bounds / Poisk optimalnyh granic vtorogo profilya
    int otstup = params.otstup;  // Search range / Diapozon poiska
    double minraz = std::numeric_limits<double>::max();
    int best_XL2 = XL1, best_XR2 = XR1;

    // Clear yP2 before each search / Ochishchaem yP2 pered kazhdym poiskom
    for(int XL2 = XL1; XL2 > XL1 - otstup; XL2--) {
        for(int XR2 = XR1; XR2 < XR1 + otstup; XR2++) {

            // Skip invalid bounds / Propuskaem nevalidnye granicy
            if(XL2 < 0 || XR2 >= n_myu || XL2 >= XR2) continue;

            // Prepare second profile / Podgotavlivaem vtoroy profil
            QVector<QPointF> yP2;
            yP2.reserve(n_myu);
            for(int i = 0; i < n_myu; i++) {
                yP2.push_back({ 1. * i, prof2[i]} );
            }

            // Reinterpolate second profile with current bounds / Pereinterpoliruem vtoroy profil s tekushchimi granicami
            QVector<QPointF> yyP2 = ImageProcessor::reInterpolateProfile(yP2, XL2, XR2, NN);
            QVector<double> yy2;
            yy2.reserve(NN);
            for(int s1 = 0; s1 < NN; s1++) {
                yy2.push_back(yyP2[s1].y());
            }

            // Normalize both profiles / Normalizuem oba profilya
            double m1 = 0.0, m2 = 0.0, D1 = 0.0, D2 = 0.0;
            for(int i = 0; i < NN; i++) {
                m1 += yy1[i];
                m2 += yy2[i];
            }
            m1 /= NN + 1;
            m2 /= NN + 1;

            for(int i = 0; i < NN; i++) {
                D1 += (yy1[i] - m1) * (yy1[i] - m1);
                D2 += (yy2[i] - m2) * (yy2[i] - m2);
            }
            D1 = sqrt(D1 / (NN + 1));
            D2 = sqrt(D2 / (NN + 1));

            // Avoid division by zero / Izbegaem deleniya na nol
            if (D1 < 1e-10) D1 = 1.0;
            if (D2 < 1e-10) D2 = 1.0;

            QVector<double> yyn1, yyn2;
            yyn1.reserve(NN);
            yyn2.reserve(NN);
            for(int s1 = 0; s1 < NN; s1++) {
                yyn1.push_back((yy1[s1] - m1) / D1);
                yyn2.push_back((yy2[s1] - m2) / D2);
            }

            // Calculate residual (sum of squared differences) / Vychislyaem nevyazku (summu kvadratov raznostey)
            double MINRAZ = 0.0;
            for(int i = 0; i < NN; i++) {
                MINRAZ += (yyn1[i] - yyn2[i]) * (yyn1[i] - yyn2[i]);
            }

            // Keep best bounds / Sohranyaem luchshie granicy
            if(minraz > MINRAZ) {
                minraz = MINRAZ;
                best_XL2 = XL2;
                best_XR2 = XR2;
            }
        }
    }

    // Final processing with optimal bounds / Finalnaya obrabotka s optimalnymi granicami
    QVector<QPointF> yP2_final;
    yP2_final.reserve(n_myu);
    for(int i = 0; i < n_myu; i++) {
        yP2_final.push_back({ 1. * i, prof2[i]} );
    }

    QVector<QPointF> yyP2_final = ImageProcessor::reInterpolateProfile(yP2_final, best_XL2, best_XR2, NN);
    QVector<double> yy2_final;
    yy2_final.reserve(NN);
    for(int s1 = 0; s1 < NN; s1++) {
        yy2_final.push_back(yyP2_final[s1].y());
    }

    // Final normalization / Finalnaya normalizaciya
    double m1_f = 0.0, m2_f = 0.0, D1_f = 0.0, D2_f = 0.0;
    for(int i = 0; i < NN; i++) {
        m1_f += yy1[i];
        m2_f += yy2_final[i];
    }
    m1_f /= NN + 1;
    m2_f /= NN + 1;

    for(int i = 0; i < NN; i++) {
        D1_f += (yy1[i] - m1_f) * (yy1[i] - m1_f);
        D2_f += (yy2_final[i] - m2_f) * (yy2_final[i] - m2_f);
    }
    D1_f = sqrt(D1_f / (NN + 1));
    D2_f = sqrt(D2_f / (NN + 1));

    // Avoid division by zero in normalization / Izbegaem deleniya na nol pri normalizacii
    if (D1_f < 1e-10) D1_f = 1.0;
    if (D2_f < 1e-10) D2_f = 1.0;

    QVector<double> yyn1_final, yyn2_final;
    yyn1_final.reserve(NN);
    yyn2_final.reserve(NN);
    for(int s1 = 0; s1 < NN; s1++) {
        yyn1_final.push_back((yy1[s1] - m1_f) / D1_f);
        yyn2_final.push_back((yy2_final[s1] - m2_f) / D2_f);
    }

    // Calculate final residual / Vychislyaem finalnuyu nevyazku
    double MINRAZ_final = 0.0;
    for(int i = 0; i < NN; i++) {
        MINRAZ_final += (yyn1_final[i] - yyn2_final[i]) * (yyn1_final[i] - yyn2_final[i]);
    }

    // Calculate shift coefficients / Vychislyaem koefficienty sdviga
    double xa = XL1, xb = XR1;
    double ya = 1. * best_XL2 - XL1;
    double yb = 1. * best_XR2 - XR1;

    // Avoid division by zero / Izbegaem deleniya na nol
    if (abs(yb - ya) < 1e-10) {
        qDebug() << "Error: Division by zero in shift calculation / Oshibka: Delenie na nol pri vychislenii sdviga";
        return result;
    }

    // Find intersection point / Nahodim tochku peresecheniya
    double x00 = (xa*yb - xb*ya) / (yb - ya);

    // Convert to mu space / Preobrazuem v prostranstvo mu
    double mu00 = -n_sigma*sigma_myu + x00*(2.*n_sigma*sigma_myu)/n_myu;

    // Final shift and scale coefficients / Finalnye koefficienty sdviga i masshtaba
    double FFF1_1 = mu00;  // Shift coefficient / Koefficient sdviga
    double FFF1_0 = (yb - ya) / (xb - xa);  // Scale coefficient / Koefficient masshtaba

    // Calculate new refined position / Vychislyaem novuyu utochnennuyu poziciyu
    double n_new = n0 + params.ex * FFF1_1;
    double m_new = m0 + params.ey * FFF1_1;

    // Fill result structure / Zapolnyaem strukturu rezultata
    result.FFF1_0 = FFF1_0;
    result.FFF1_1 = FFF1_1;
    result.refinedPosition = QPointF(n_new, m_new);
    result.residual = MINRAZ_final;
    result.success = true;

    return result;
}
