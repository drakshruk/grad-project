#include "imageprocessor.h"

/*
// * EN:
// *  substracts Image2 from Image1
// * RU:
// *  vichitaet Image2 iz Image1
// *
//const QImage ImageProcessor::substractTwoImages(const QImage &Im1, const QImage &Im2)
//{
//    QImage res(Im1.width(),Im1.height(),QImage::Format_ARGB32);
//    for(int i = 0; i < Im1.width(); i++){
//        for(int j = 0; j < Im1.height(); j++){
//            res.setPixel(i,j,qRgb(abs(qRed(Im1.pixel(i,j)) - qRed(Im2.pixel(i,j))),
//                                   abs(qGreen(Im1.pixel(i,j)) - qGreen(Im2.pixel(i,j))),
//                                   abs(qBlue(Im1.pixel(i,j)) - qBlue(Im2.pixel(i,j)))));
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  substracts ddMat2 from ddMat1 with 0 sampling
// * RU:
// *  vichitaet ddMat2 iz ddMat1
// *
//const Matrix2D<double> ImageProcessor::substractTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
////            if(res[i][j] < 0.0) res[i][j] = 0.0; //
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  adds ddMat2 to ddMat1 with 255 sampling
// * RU:
// *  pribavlyaet ddMat2 k ddMat1
// *
//const Matrix2D<double> ImageProcessor::addTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            res[i][j] = ddMat1[i][j] + ddMat2[i][j];
//            if(res[i][j] > 255.0) res[i][j] = 255.0; //
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  multiples ddMat1 by ddMat2 with 255 sampling
// * RU:
// *
// *
//const Matrix2D<double> ImageProcessor::multiplyTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            res[i][j] = ddMat1[i][j] * ddMat2[i][j];
//            if(res[i][j] > 255.0) res[i][j] = 255.0; //
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  divides ddMat1 by ddMat2
// * RU:
// *
// *
//const Matrix2D<double> ImageProcessor::divideTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            if(abs(ddMat2[i][j]) < 1e-6) res[i][j] = 255.;
//            else res[i][j] = ddMat1[i][j] / ddMat2[i][j];
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  sets pixels to min between ddMat1 and ddMat2
// * RU:
// *
// *
//const Matrix2D<double> ImageProcessor::minTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            res[i][j] = min(ddMat1[i][j],ddMat2[i][j]);
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  sets pixels to max between ddMat1 and ddMat2
// * RU:
// *
// *
//const Matrix2D<double> ImageProcessor::maxTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            res[i][j] = max(ddMat1[i][j],ddMat2[i][j]);
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  sets pixels to ddMat1 AND ddMat2
// * RU:
// *
// *
//const Matrix2D<double> ImageProcessor::AND_TwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            res[i][j] = static_cast<int>(ddMat1[i][j]) & static_cast<int>(ddMat2[i][j]);
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  sets pixels to ddMat1 OR ddMat2
// * RU:
// *
// *
//const Matrix2D<double> ImageProcessor::OR_TwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            res[i][j] = static_cast<int>(ddMat1[i][j]) | static_cast<int>(ddMat2[i][j]);
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  sets pixels to ddMat1 XOR ddMat2
// * RU:
// *
// *
//const Matrix2D<double> ImageProcessor::XOR_TwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            res[i][j] = static_cast<int>(ddMat1[i][j]) ^ static_cast<int>(ddMat2[i][j]);
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  sets pixels to average between ddMat1 and ddMat2
// * RU:
// *
// *
//const Matrix2D<double> ImageProcessor::averageTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            res[i][j] = (ddMat1[i][j] + ddMat2[i][j]) / 2.;
//        }
//    }
//    return res;
//}

///*
// * EN:
// *  sets pixels to absolute difference between ddMat1 and ddMat2
// * RU:
// *
// *
//const Matrix2D<double> ImageProcessor::differenceTwoMatrices(const Matrix2D<double> &ddMat1, const Matrix2D<double> &ddMat2)
//{
//    size_t w = min(ddMat1.size(), ddMat2.size());
//    size_t h = min(ddMat1[0].size(), ddMat2[0].size());
//    Matrix2D<double> res(ddMat1.size(), std::vector<double>(ddMat1[0].size(), 0.0));
//    for(size_t i = 0; i < w; i++){
//        for(size_t j = 0; j < h; j++){
//            res[i][j] = ddMat1[i][j] - ddMat2[i][j];
//            res[i][j] = abs(ddMat1[i][j] - ddMat2[i][j]);
//            if(res[i][j] < 0) res[i][j] = -res[i][j];
//        }
//    }
//    return res;
//}
*/

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

        if (xCoord < 0 || xCoord >= ddImMat.size() || yCoord < 0 || yCoord >= ddImMat[0].size()) {
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

    if (convCore.size() != kernelSize || convCore[0].size() != kernelSize) {
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
 *
 * RU:
 *
 */
const std::pair<int,int> ImageProcessor::findProfileBounds(const QVector<QPointF> profile, int step)
{
    if (profile.size() < 3) {
            return {0, profile.size() - 1};
    }

    int center = profile.size() / 2;
    int leftBound = center, rightBound = center;
    bool leftIncr = true, rightIncr = true, decrLeft = (profile[leftBound-1].y() < profile[leftBound].y());
    while(leftIncr || rightIncr){
        if(leftBound <= 0 || rightBound >= profile.size()-1) break;
        if(decrLeft){
            if(profile[leftBound-1].y() < profile[leftBound].y() && leftIncr){
                leftBound--;
            }
            else if(leftIncr){
                if(leftBound < step) leftBound = 0;
                else leftBound -= step;
                leftIncr = false;
            }
            if(profile[rightBound+1].y() > profile[rightBound].y() && rightIncr){
                rightBound++;
            }
            else if(rightIncr){
                if(rightBound > profile.size() - step - 1) rightBound = profile.size() - 1;
                else rightBound += step;
                rightIncr = false;
            }
        }
        else{
            if(profile[leftBound-1].y() > profile[leftBound].y() && leftIncr){
                leftBound--;
            }
            else if(leftIncr){
                if(leftBound < step) leftBound = 0;
                else leftBound -= step;
                leftIncr = false;
            }
        }
        if(profile[rightBound+1].y() < profile[rightBound].y() && rightIncr){
            rightBound++;
        }
        else if(rightIncr){
            if(rightBound > profile.size() - step - 1) rightBound = profile.size() - 1;
            else rightBound += step;
            rightIncr = false;
        }
    }
    return {leftBound, rightBound};
}

/*
 * EN:
 *
 * RU:
 *
 */
std::pair<int,int> ImageProcessor::alignProfiles(QVector<QPointF> &profileForReference, QVector<QPointF> &profileToAlign, int left, int right, int searchRange)
{
    if(profileForReference.empty() || profileToAlign.empty()) return {-1,-1};

    std::pair<int,int> optimalShift = {0,0};

    double minResidual = std::numeric_limits<double>::max();

    double median1 = 0.0;
    double dispersion1 = 0.0;

    for(int i = 0; i < profileForReference.size(); i++){
        median1 += profileForReference[i].y();
    }
    median1 /= profileForReference.size();

    for(int i = 0; i < profileForReference.size(); i++){
        dispersion1 += (profileForReference[i].y()-median1)*(profileForReference[i].y()-median1);
    }
    dispersion1 = sqrt(dispersion1)/profileForReference.size();

    for(int i = 0; i < profileForReference.size(); i++){
        profileForReference[i].setY((profileForReference[i].y()-median1)/dispersion1);
    }
    int leftEdge = left - searchRange;
    leftEdge = leftEdge < 0 ? 0 : leftEdge;

    int rightEdge = right + searchRange;
    rightEdge = rightEdge >= profileToAlign.size() ? profileToAlign.size() - 1 : rightEdge;

    for(int leftIncr = leftEdge; leftIncr <= left; leftIncr++){
        for(int rightIncr = rightEdge; rightIncr >= right; rightIncr--){

            QVector<QPointF> shiftedProfile = cutProfile(profileToAlign, leftIncr, rightIncr);
            shiftedProfile = reInterpolateProfile(shiftedProfile, profileForReference.size());

            double median2 = 0.0;
            double dispersion2 = 0.0;

            for(int i = 0; i < profileToAlign.size(); i++){
                median2 += profileToAlign[i].y();
            }
            median2 /= profileToAlign.size();

            for(int i = 0; i < profileToAlign.size(); i++){
                dispersion2 += (profileToAlign[i].y()-median2)*(profileToAlign[i].y()-median2);
            }
            dispersion2 = sqrt(dispersion2)/profileToAlign.size();

            for(int i = 0; i < profileToAlign.size(); i++){
                profileToAlign[i].setY((profileToAlign[i].y()-median2)/dispersion2);
            }

            double tmp = calculateResidualOfProfiles(profileForReference,shiftedProfile);
            if(tmp < minResidual) {
                optimalShift = {leftIncr,rightIncr};
                minResidual = tmp;
            }
        }
    }
    return optimalShift;
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
 *  reinterpolates profile to the new length 
 *  using linear interpolation
 * RU:
 *  pereinterpolirayet profil' do novoy dliny
 *  ispol'zuya lineynuyu interpolyatsiyu
 */
const QVector<QPointF> ImageProcessor::reInterpolateProfile(const QVector<QPointF> &profile, int newLen)
{
    QVector<QPointF> newProfile;

    newProfile.clear();
    newProfile.reserve(newLen);

    if (profile.size() < 2 || newLen <= 0) {
        return profile;
    }

    double xL1 = profile.front().x();
    double xR1 = profile.back().x();

    newProfile.reserve(newLen);

    double dxNew = (xR1 - xL1) / (newLen - 1);

    for (int i = 0; i < newLen; i++) {
        double x = xL1 + i * dxNew;

        int segmentIndex = 0;
        while (segmentIndex < profile.size() - 1 && profile[segmentIndex + 1].x() < x) {
            segmentIndex++;
        }

        if (segmentIndex >= profile.size() - 1) {
            newProfile.append(QPointF(profile.size()-1,profile.last().y()));
            continue;
        }

        const QPointF& p1 = profile.at(segmentIndex);
        const QPointF& p2 = profile.at(segmentIndex + 1);

        if (p2.x() == p1.x()) {
            newProfile.append(QPointF(i, p1.y()));
        } else {
            double y = p1.y() + (x - p1.x()) * (p2.y() - p1.y()) / (p2.x() - p1.x());
            newProfile.append(QPointF(i, y));
        }
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
    int rows = ddImageMat.size(), cols = ddImageMat[0].size();
    int iRad = ddConvCore.size();
    Matrix2D<double> ddResMat(rows, std::vector<double>(cols));
    double matSum = 0., convSum = 0.;
    for(int ii = 0; ii < iRad; ii++) {
        for(int jj = 0; jj < iRad; jj++) {
            convSum += ddConvCore[ii][jj];
        }
    }
   if (convSum <= 0.) convSum = 1.0;

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            matSum = 0.;
            for(int ii = 0; ii < iRad; ii++)
            {
                int it = abs(i - iRad/2 + ii);
                if( it >= rows)
                {
                    it = 2*(rows - 1) - it;
                }
                for(int jj = 0; jj < iRad; jj++)
                {
                    int jt = abs(j - iRad/2 + jj);

                    if( jt >= cols)
                    {
                        jt = 2*(cols - 1) - jt;
                    }
                    matSum += ddImageMat[it][jt] * ddConvCore[ii][jj];
                }
            }
            matSum /= convSum;
            ddResMat[i][j] = matSum;
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
    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize));
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    double sum = 0.;
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = abs(xCenter - i);
            int dy = abs(yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = 1.0 / (dSigma * sqrt(2.0 * pi)) * exp((-1.0 / 2.0) * (Rad * Rad) / (dSigma * dSigma));
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
    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize));
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    double sum = 0.0;
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = (xCenter - i);
            int dy = (yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = dx / (dSigma * sqrt(2.0 * pi)) * exp((-1.0 / 2.0) * (Rad * Rad) / (dSigma * dSigma));
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
    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize));
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    double sum = 0.0;
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = (xCenter - i);
            int dy = (yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = dy / (dSigma * sqrt(2.0 * pi)) * exp((-1.0 / 2.0) * (Rad * Rad) / (dSigma * dSigma));
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
    Matrix2D<double> ddRes(iXsize, std::vector<double>(iYsize));
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    double sum = 0.0;
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = abs(xCenter - i);
            int dy = abs(yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = (Rad * Rad / dSigma / dSigma - 2.) * exp( -1. * Rad * Rad / 2 / dSigma / dSigma);
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
    Matrix2D<double> mat1;
    Matrix2D<double> mat2;

    mat2 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma));
    mat1 = convMat(ddImageMat, getGauss(iRad, iRad, dSigma*2.0));
    mat1 = elementWiseOperation(mat1, mat2, MatrixLambdas::Subtract<double>{});
    ddMatForProfile = mat1;
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
