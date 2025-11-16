#include "../include/mainwindow.h"


// math functions

double min(double** ddMat, int iXsize, int iYsize);
double max(double** ddMat, int iXsize, int iYsize);
int min(int** ddMat, int iXsize, int iYsize);
int max(int** ddMat, int iXsize, int iYsize);
QImage substractIm(QImage Im1, QImage Im2);
bool containsVP(QVector<QPoint>* vec, QPoint p);
int** selectEdge(QPoint pos, int** attMat, QVector<QPoint> *selectedEdge);
double** substractMat(double** ddMat1, double** ddMat2, int iXsize, int iYsize);

template<typename T>
void freeMatrix(T** matrix, int rows) {
    if (matrix == nullptr) return;

    for (int i = 0; i < rows; i++) {
        if (matrix[i] != nullptr) {
            delete[] matrix[i];
            matrix[i] = nullptr;
        }
    }
    delete[] matrix;
}

template<typename T>
void copyMatrix(T** src, T** dest, int rows, int cols)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

double sum(QVector<QPointF> vec){
    double res = 0.;
    for(int i = 0; i < vec.size(); i++){
        res += vec[i].y();
    }
    return res;
}

// profile functions
QVector<QPointF> buildProfile1(int x0, int y0, double** ddImMat, double xGrad, double yGrad, int xSize, int ySize, double dSigma);
QVector<QPointF> buildProfile2(int x0, int y0, double** ddImMat, double xGrad, double yGrad, int xSize, int ySize, double dSigma, int nSig, QVector<QPoint>* indices);
QVector<QPointF> cutProfile(QVector<QPointF> profile); // cuts profile from closest to middle local min to closest to middle local max
QVector<QPointF> cutProfile(QVector<QPointF> profile, int left, int right);
void alignProfiles(QVector<QPointF> profile1, QVector<QPointF> profile2, int len);
double differenceOfTwoProfiles(QVector<QPointF> profile1, QVector<QPointF> profile2);
void reInterpolateProfile(QVector<QPointF> *profile, int newLen, QVector<QPointF> *newProfile);

// to/from/combine image functions
double** fromGrayImage(QImage image);
QImage toGrayImage(double** ddMat, int iXsize, int iYsize);
QImage toBlueRedImage(double** ddMat, int iXsize, int iYsize, double dRedMax, double dBlueMax);
QImage combineWithEdge(double** ddMat, double** ddEdge, int iXsize, int iYsize);
QImage combineImageWithProfile(QImage im1, QVector<QPoint> profile);

//convolutional functions
QImage convImage(QImage image, double** ddLapl, int iRad);
double** convMat(double** ddMat, double** ddConvCore, int iRad, int iXsize, int iYsize);
double** getGauss(int iXsize, int iYsize, double dSigma);
double** getXGradCore(int iXsize, int iYsize, double dSigma);
double** getYGradCore(int iXsize, int iYsize, double dSigma);
double** getLapl(int iXsize, int iYsize, double dSigma);
double** findEdges(double** ddMat, int iXsize, int iYsize, int iRad, int** attMat);
double** findEdgesRB(double** ddMat, int iXsize, int iYsize, int iRad);
QImage gaussianEdgeDetection(double** ddMat, int iXsize, int iYsize, double dSigma, int iRad, double** matForProfile);


//samples
const QImage sampleTwoHollows()
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

const QImage sampleTwoHollowsBig()
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
//

double min(double** ddMat, int iXsize, int iYsize)
{
    double min = ddMat[0][0];
    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            if(min > ddMat[i][j]) min = ddMat[i][j];
        }
    }
    return min;
}

double max(double** ddMat, int iXsize, int iYsize)
{
    double max = ddMat[0][0];
    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            if(max < ddMat[i][j]) max = ddMat[i][j];
        }
    }
    return max;
}

int min(int** ddMat, int iXsize, int iYsize)
{
    int min = ddMat[0][0];
    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            if(min > ddMat[i][j]) min = ddMat[i][j];
        }
    }
    return min;
}

int max(int** ddMat, int iXsize, int iYsize)
{
    int max = ddMat[0][0];
    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            if(max < ddMat[i][j]) max = ddMat[i][j];
        }
    }
    return max;
}

QImage substractIm(QImage Im1, QImage Im2)
{
    for(int i = 0; i < Im1.width(); i++)
    {
        for(int j = 0; j < Im1.height(); j++)
        {
            Im1.setPixel(i,j,qRgb(abs(qRed(Im1.pixel(i,j)) - qRed(Im2.pixel(i,j))),
                                  abs(qGreen(Im1.pixel(i,j)) - qGreen(Im2.pixel(i,j))),
                                  abs(qBlue(Im1.pixel(i,j)) - qBlue(Im2.pixel(i,j)))));
        }
    }
    return Im1;
}

QImage toGrayImage(double** ddMat, int iXsize, int iYsize)
{
    QImage resIm(iXsize, iYsize, QImage::Format_ARGB32);

    double dMin = min(ddMat, iXsize, iYsize);
    double dMax = max(ddMat, iXsize, iYsize);

    double** ddTmp = new double*[iXsize];
    for(int i = 0; i < iXsize; i++)
    {
        ddTmp[i] = new double[iYsize];
        for(int j = 0; j < iYsize; j++) ddTmp[i][j] = ddMat[i][j];
    }

    if(dMin < 0)
    {
        for(int i = 0; i < iXsize; i++)
        {
            for(int j = 0; j < iYsize; j++)
            {
                ddTmp[i][j] += dMin;
                ddTmp[i][j] *= (255 / (dMax + dMin));
                resIm.setPixel(i,j,qRgb((int)ddTmp[i][j], (int)ddTmp[i][j],(int)ddTmp[i][j]));
            }
        }
    }
    if(dMin >= 0)
    {
        for(int i = 0; i < iXsize; i++)
        {
            for(int j = 0; j < iYsize; j++)
            {
                ddTmp[i][j] *= (255 / dMax);
                resIm.setPixel(i,j,qRgb((int)ddTmp[i][j], (int)ddTmp[i][j],(int)ddTmp[i][j]));
            }
        }
    }
    return resIm;
}

QImage toBlueRedImage(double** ddMat, int iXsize, int iYsize, double dRedMax, double dBlueMax)
{
    QImage resIm(iXsize, iYsize, QImage::Format_ARGB32);

    double dMin = min(ddMat, iXsize, iYsize);
    double dMax = max(ddMat, iXsize, iYsize);

    double** ddTmp = new double*[iXsize];
    for(int i = 0; i < iXsize; i++)
    {
        ddTmp[i] = new double[iYsize];
        for(int j = 0; j < iYsize; j++) ddTmp[i][j] = ddMat[i][j];
    }

    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            if(ddMat[i][j] < 0)
            {
                ddTmp[i][j] *= (dBlueMax/dMin);
                resIm.setPixel(i,j,qRgb(0, 0,(int)ddTmp[i][j]));
            }
            else if(ddMat[i][j] > 0)
            {
                ddTmp[i][j] *= (dRedMax/dMax);
                resIm.setPixel(i,j,qRgb((int)ddTmp[i][j], 0, 0));
            }
            else resIm.setPixel(i,j,qRgb(0,0,0));
        }
    }
    return resIm;
}

QImage combineWithEdge(double** ddMat, double** ddEdge, int iXsize, int iYsize)
{
    QImage resIm(iXsize, iYsize, QImage::Format_ARGB32);
    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            if(ddEdge[i][j] == 255) resIm.setPixel(i,j,qRgb(255,255,255));
            else resIm.setPixel(i,j,qRgb((int)ddMat[i][j], (int)ddMat[i][j], (int)ddMat[i][j]));
        }
    }
    return resIm;
}

QImage combineImageWithProfile(QImage im1, QVector<QPoint> profile)
{
    QImage resIm(im1.width(), im1.height(),QImage::Format_ARGB32);
    for(int i = 0; i < profile.size(); i++){
        im1.setPixel(profile[i].x(),profile[i].y(),qRgb(0,255,255));
    }
    return resIm;
}

double** fromGrayImage(QImage image)
{
    double** ddRes = new double*[image.width()];
    for(int i = 0; i < image.width(); i++) ddRes[i] = new double[image.height()];
    for(int i = 0; i < image.width(); i++)
    {
        for(int j = 0; j < image.height(); j++)
        {
            ddRes[i][j] = 1.0 * qGray(image.pixel(i,j));
        }
    }
    return ddRes;
}

double** substractMat(double** ddMat1, double** ddMat2, int iXsize, int iYsize)
{
    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            ddMat1[i][j] -= ddMat2[i][j];
        }
    }
    return ddMat1;
}

bool containsVP(QVector<QPoint>* vec, QPoint p)
{
    for(int i = 0; i < vec->size(); i++)
    {
        if(vec->at(i).x() == p.x() && vec->at(i).y() == p.y()) return true;
    }
    return false;
}

int** selectEdge(QPoint pos, int** attMat, QVector<QPoint> *selectedEdge)
{
    if(!selectedEdge->empty()) selectedEdge->clear();
    if(attMat[pos.x()][pos.y()] == int(attribute::isEdge))
    {
        selectedEdge->push_back(pos);
        attMat[pos.x()][pos.y()] = int(attribute::isSelectedEdge);
        size_t prevSize = 0, newSize = 1;
        while(prevSize < newSize)
        {
            prevSize = selectedEdge->size();
            for(int s = 0; s < selectedEdge->size(); s++)
            {
                for(int ii = -1; ii <= 1; ii++)
                {
                    for(int jj = -1; jj <= 1; jj++)
                    {
                        if(attMat[selectedEdge->at(s).x() + ii][selectedEdge->at(s).y() + jj] == int(attribute::isEdge) &&
                                !containsVP(selectedEdge,QPoint(selectedEdge->at(s).x() + ii, selectedEdge->at(s).y() + jj)))
                        {
                            selectedEdge->push_back(QPoint(selectedEdge->at(s).x() + ii, selectedEdge->at(s).y() + jj));
                            attMat[selectedEdge->at(s).x() + ii][selectedEdge->at(s).y() + jj] = int(attribute::isSelectedEdge);
                        }
                    }
                }
            }
            newSize = selectedEdge->size();
        }
    }
    return attMat;
}

QVector<QPointF> buildProfile1(int x0, int y0, double** ddImMat, double xGrad, double yGrad, int xSize, int ySize, double dSigma)
{
    QVector<QPointF> res;
    double nsig = 3, n = int(nsig*dSigma);
    double muS, x, y;
//    qDebug() << n;
    for(int s = 0; s < n; s++){
//        qDebug() << s;
        muS = -nsig*dSigma + 2.0*nsig*dSigma*s/n;
        x = x0 + muS*xGrad;
        y = y0 + muS*yGrad;
//        qDebug() << int(x+0.5) << int(y+0.5);
//        qDebug() << ddImMat[int(x+0.5)][int(y+0.5)];
        if(x < 0 || x >= xSize || y < 0 || y >= ySize) res.push_back(QPoint(0.,0.));
        else res.push_back(QPointF(s,ddImMat[int(x+0.5)][int(y+0.5)]));
    }
    return res;
}

QVector<QPointF> buildProfile2(int x0, int y0, double** ddImMat, double xGrad, double yGrad, int xSize, int ySize, double dSigma, int nSig, QVector<QPoint>* indices)
{
    QVector<QPointF> res;
    int iRad = 2;
    double n = int(nSig*dSigma);
    double muS, x, y, yProf;
    double** ddLapl = getLapl(5,5,dSigma);
//    qDebug() << n;
    for(int s = 0; s < n; s++){
//        qDebug() << s;
        muS = -nSig*dSigma + 2.0*nSig*dSigma/n*s;
        x = x0 + muS*xGrad;
        y = y0 + muS*yGrad;
//        qDebug() << int(x+0.5) << int(y+0.5);
//        qDebug() << ddImMat[int(x+0.5)][int(y+0.5)];
        if(x < 0 || x >= xSize || y < 0 || y >= ySize) res.push_back(QPoint(0.,0.));
        else
        {
            yProf = 0.;
            for(int i = -iRad; i <= iRad; i++){
                for(int j = -iRad; j < iRad; j++){
                    yProf += ddLapl[i+iRad][j+iRad]*ddImMat[int(x+0.5+i)][int(y+0.5+j)];
                }
            }
            indices->push_back(QPoint(int(x+0.5),int(y+0.5)));
            res.push_back(QPointF(s,ddImMat[int(x+0.5)][int(y+0.5)]));
        }
    }
    return res;
}

QVector<QPointF> cutProfile(QVector<QPointF>* profile, int otstup) //  goes from the center to left/ right, finds min/max and goes another 5~ pixels
{
    QVector<QPointF> res;
    int left = profile->size()/2, right = profile->size()/2;
    bool leftIncr = true, rightIncr = true, decrLeft = (profile->at(left-1).y() < profile->at(left).y());
    qDebug() << decrLeft << "123123";
    while(leftIncr || rightIncr){
        if(left <= 0 || right >= profile->size()-1) break;
        if(decrLeft){
            if(profile->at(left-1).y() < profile->at(left).y() && leftIncr){
                left--;
            }
            else if(leftIncr){
                if(left < otstup) left = 0;
                else left -= otstup;
                leftIncr = false;
            }
            if(profile->at(right+1).y() > profile->at(right).y() && rightIncr){
                right++;
            }
            else if(rightIncr){
                if(right > profile->size() - otstup - 1) right = profile->size() - 1;
                else right += otstup;
                rightIncr = false;
            }
        }
        else{
            if(profile->at(left-1).y() > profile->at(left).y() && leftIncr){
                left--;
                qDebug()  << "left: " << left << " " << profile->at(left).y();
                qDebug() << leftIncr;
            }
            else if(leftIncr){
                if(left < otstup) left = 0;
                else left -= otstup;
                qDebug()  << "left: " << left << " " << profile->at(left).y();
                qDebug() << leftIncr;
                leftIncr = false;
            }
        }
        if(profile->at(right+1).y() < profile->at(right).y() && rightIncr){
            right++;
            qDebug()  << "right: " << right << " " << profile->at(right).y();
            qDebug() << rightIncr;
        }
        else if(rightIncr){
            if(right > profile->size() - otstup - 1) right = profile->size() - 1;
            else right += otstup;
            qDebug()  << "right: " << right << " " << profile->at(right).y();
            qDebug() << rightIncr;
            rightIncr = false;
        }
    }

    for(int i = left; i < right; i++){
        res.push_back(profile->at(i));
    }
    qDebug() << decrLeft << " " << left << " " << right;
    return res;
}

QVector<QPointF> cutProfile(QVector<QPointF> profile, int left, int right)
{
    QVector<QPointF> res;
    if(right >= profile.size()) right = profile.size()-1;
    if(left < 0) left = 0;
    for(int i = left; i <= right; i++){
        res.append(profile[i]);
    }
    return res;
}

void alignProfiles(QVector<QPointF> profile1, QVector<QPointF> profile2, int len, double dSigma1, double dSigma2)
{
    QVector<QPointF> reIntProfile1, reIntProfile2, cutProfile2;
    int xL1 = profile1.first().x(), xR1 = profile1.last().x();
    int xL2, xR2, minXL2 = 0, minXR2 = 0;
    double delta = INFINITY, shift = 20;
    double m1, m2, D1, D2;
    profile1 = cutProfile(&profile1, len);
    xL2 = xL1 - 1;
    xR2 = xR1 + 1;
    reInterpolateProfile(&profile1, len, &reIntProfile1);
    for(int i = 0; i < shift; i++){
        for(int j = 0; j < shift; j++){
            cutProfile2.clear();
            cutProfile2 = cutProfile(profile2, xL2, xR2);
            reInterpolateProfile(&cutProfile2, len, &reIntProfile2);
            m1 = 1./(len+1.) * sum(reIntProfile1);
            m2 = 1./(len+1.) * sum(reIntProfile2);

            double tmp = differenceOfTwoProfiles(reIntProfile1, reIntProfile2);
            if(delta > tmp){
                delta = tmp;
                minXL2 = xL2;
                minXR2 = xR2;
            }
        }
    }


}

double differenceOfTwoProfiles(QVector<QPointF> profile1, QVector<QPointF> profile2)
{
    double res = 0;
    for(int i = 0; i < profile1.size(); i++){
        res += (profile1[i].y()-profile2[i].y())*(profile1[i].y()-profile2[i].y());
    }
    return res;
}

void reInterpolateProfile(QVector<QPointF> *profile, int newLen, QVector<QPointF> *newProfile)
{
    if (profile->size() < 2 || newLen <= 0) {
        newProfile->clear();
        return;
    }

    double xL1 = profile->front().x();
    double xR1 = profile->back().x();

    newProfile->clear();
    newProfile->reserve(newLen);

    double dxNew = (xR1 - xL1) / (newLen - 1);

    for (int i = 0; i < newLen; i++) {
        double x = xL1 + i * dxNew;

        int segmentIndex = 0;
        while (segmentIndex < profile->size() - 1 && profile->at(segmentIndex + 1).x() < x) {
            segmentIndex++;
        }

        if (segmentIndex >= profile->size() - 1) {
            newProfile->append(QPointF(profile->size()-1,profile->last().y()));
            continue;
        }

        const QPointF& p1 = profile->at(segmentIndex);
        const QPointF& p2 = profile->at(segmentIndex + 1);

        if (p2.x() == p1.x()) {
            newProfile->append(QPointF(i, p1.y()));
        } else {
            double y = p1.y() + (x - p1.x()) * (p2.y() - p1.y()) / (p2.x() - p1.x());
            newProfile->append(QPointF(i, y));
        }
    }
}

QImage convImage(QImage image, double** ddLapl, int iRad)
{
    QImage resIm = image;
    for(int i = 0; i < image.width(); i++)
    {
        for(int j = 0; j < image.height(); j++)
        {
            double rSum = 0, gSum = 0, bSum = 0, lSum = 0;
            for(int ii = 0; ii < iRad; ii++)
            {
                for(int jj = 0; jj < iRad; jj++)
                {
                    int it = abs(i - iRad/2 + ii);
                    int jt = abs(j - iRad/2 + jj);
                    if( it >= image.width())
                    {
                        it = 2*(image.width() - 1) - it;
                    }
                    if( jt >= image.height())
                    {
                        jt = 2*(image.height() - 1) - jt;
                    }
                    rSum += (double)qRed(image.pixel(it,jt)) * ddLapl[ii][jj];
                    gSum += (double)qGreen(image.pixel(it,jt)) * ddLapl[ii][jj];
                    bSum += (double)qBlue(image.pixel(it,jt)) * ddLapl[ii][jj];
                    lSum += ddLapl[ii][jj];
                }
            }
            if (lSum <= 0) lSum = 1;
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

double** convMat(double** ddMat, double** ddConvCore, int iRad, int iXsize, int iYsize)
{
    double** ddRes = new double*[iXsize];
    double matSum = 0, convSum = 0;
    for(int i = 0; i < iXsize; i++) ddRes[i] = new double[iYsize];

    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            matSum = 0; convSum = 0;
            for(int ii = 0; ii < iRad; ii++)
            {
                for(int jj = 0; jj < iRad; jj++)
                {
                    int it = abs(i - iRad/2 + ii);
                    int jt = abs(j - iRad/2 + jj);
                    if( it >= iXsize)
                    {
                        it = 2*(iXsize - 1) - it;
                    }
                    if( jt >= iYsize)
                    {
                        jt = 2*(iYsize - 1) - jt;
                    }
                    matSum += ddMat[it][jt] * ddConvCore[ii][jj];
                    convSum += ddConvCore[ii][jj];
                }
            }
            if (convSum <= 0) convSum = 1;
            matSum /= convSum;
            ddRes[i][j] = matSum;
        }
    }
    return ddRes;
}

double** getGauss(int iXsize, int iYsize, double dSigma)
{
    double** ddRes = new double*[iXsize];
    for(int i = 0; i < iXsize; i++) ddRes[i] = new double[iYsize];
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = abs(xCenter - i);
            int dy = abs(yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = 1.0 / (dSigma * sqrt(2.0 * pi)) * exp((-1.0 / 2.0) * (Rad * Rad) / (dSigma * dSigma));
        }
    }
    return ddRes;
}

double** getXGradCore(int iXsize, int iYsize, double dSigma)
{
    double** ddRes = new double*[iXsize];
    for(int i = 0; i < iXsize; i++) ddRes[i] = new double[iYsize];

    int xCenter = iXsize / 2; int yCenter = iYsize / 2;

    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = (xCenter - i);
            int dy = (yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = dx / (dSigma * sqrt(2.0 * pi)) * exp((-1.0 / 2.0) * (Rad * Rad) / (dSigma * dSigma));

        }
    }
    return ddRes;
}

double** getYGradCore(int iXsize, int iYsize, double dSigma)
{
    double** ddRes = new double*[iXsize];
    for(int i = 0; i < iXsize; i++) ddRes[i] = new double[iYsize];

    int xCenter = iXsize / 2; int yCenter = iYsize / 2;

    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = (xCenter - i);
            int dy = (yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = dy / (dSigma * sqrt(2.0 * pi)) * exp((-1.0 / 2.0) * (Rad * Rad) / (dSigma * dSigma));

        }
    }
    return ddRes;
}

double** getLapl(int iXsize, int iYsize, double dSigma)
{
    double** ddRes = new double*[iXsize];
    for(int i = 0; i < iXsize; i++) ddRes[i] = new double[iYsize];
    int xCenter = iXsize / 2; int yCenter = iYsize / 2;
    for (int i = 0; i < iXsize; i++)
    {
        for (int j = 0; j < iYsize; j++)
        {
            int dx = abs(xCenter - i);
            int dy = abs(yCenter - j);
            double Rad = sqrt(dx * dx + dy * dy);
            ddRes[i][j] = (Rad * Rad / dSigma / dSigma - 2.) * exp( -1. * Rad * Rad / 2 / dSigma / dSigma);
        }
    }
    return ddRes;
}

double** findEdges(double** ddMat, int iXsize, int iYsize, int iRad, int** attMat)
{
    double** ddEdge = new double*[iXsize];
    for(int i = 0; i < iXsize; i++) ddEdge[i] = new double[iYsize];

    bool bPos, bNeg;
    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            bPos = false; bNeg = false;
            for(int ii = 0; ii < iRad; ii++)
            {
                for(int jj = 0; jj < iRad; jj++)
                {
                    int it = abs(i - iRad/2 + ii);
                    int jt = abs(j - iRad/2 + jj);
                    if( it >= iXsize)
                    {
                        it = 2*(iXsize - 1) - it;
                    }
                    if( jt >= iYsize)
                    {
                        jt = 2*(iYsize - 1) - jt;
                    }
                    if(ddMat[it][jt] > 0) bPos = true;
                    else if(ddMat[it][jt] < 0) bNeg = true;
                }
            }
            if(bPos && bNeg)
            {
                attMat[i][j] = int(attribute::isEdge);
                ddEdge[i][j] = 255.;
            }
            else ddEdge[i][j] = 0.;
        }
    }
    return ddEdge;
}

double** findEdgesRB(double** ddMat, int iXsize, int iYsize, int iRad)
{
    double** ddEdge = new double*[iXsize];
    for(int i = 0; i < iXsize; i++) ddEdge[i] = new double[iYsize];

    bool bPos, bNeg;
    for(int i = 0; i < iXsize; i++)
    {
        for(int j = 0; j < iYsize; j++)
        {
            bPos = false; bNeg = false;
            for(int ii = 0; ii < iRad; ii++)
            {
                for(int jj = 0; jj < iRad; jj++)
                {
                    int it = abs(i - iRad/2 + ii);
                    int jt = abs(j - iRad/2 + jj);
                    if( it >= iXsize)
                    {
                        it = 2*(iXsize - 1) - it;
                    }
                    if( jt >= iYsize)
                    {
                        jt = 2*(iYsize - 1) - jt;
                    }
                    if(ddMat[it][jt] > 0) bPos = true;
                    else if(ddMat[it][jt] < 0) bNeg = true;
                }
            }
            if(bPos && bNeg)
            {
                if(ddMat[i][j] > 0) ddEdge[i][j] = 1;
                if(ddMat[i][j] <= 0) ddEdge[i][j] = -1;
            }
            else ddEdge[i][j] = 0;
        }
    }
    return ddEdge;
}

QImage gaussianEdgeDetection(double** ddMat, int iXsize, int iYsize, double dSigma, int iRad, double** matForProfile, int** attMat)
{
    double** mat1 = nullptr;
    double** mat2 = nullptr;

    mat2 = convMat(ddMat, getGauss(iRad, iRad, dSigma), iRad, iXsize, iYsize);
    mat1 = convMat(ddMat, getGauss(iRad, iRad, dSigma*2.0), iRad, iXsize, iYsize);
    mat1 = substractMat(mat1, mat2, iXsize, iYsize);
    freeMatrix(matForProfile,iXsize);
    matForProfile = new double*[iXsize];
    for(int i = 0; i < iXsize; i++){
        matForProfile[i] = new double[iYsize];
    }
    copyMatrix(mat1, matForProfile, iXsize, iYsize);
    mat1 = findEdges(mat1, iXsize, iYsize, 3, attMat);
    QImage image =  toGrayImage(mat1, iXsize, iYsize);
    freeMatrix(mat1, iXsize);
    freeMatrix(mat2, iXsize);
    return image;
}
