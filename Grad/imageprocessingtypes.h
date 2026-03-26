#ifndef IMAGEPROCESSINGTYPES_H
#define IMAGEPROCESSINGTYPES_H

#include <QImage>
#include <QVector>
#include <QPoint>
#include <cmath>
#include <QDebug>
#include "matrixlambdas.h"

template<typename T>
using Matrix3D = std::vector<std::vector<std::vector<T>>>;

template<typename T>
using Matrix2D = std::vector<std::vector<T>>;

template<typename T>
using Matrix1D = std::vector<T>;

template<typename T>
using QMatrix2D = QVector<QVector<T>>;

template<typename T>
using QMatrix1D = QVector<T>;


const double pi = 4. * atan(1);

class App_Stats {

public:
    double max_val = 0;
    double min_val = 0;
    double median_val = 0;
    double dispersion = 0;
    double sko = 0;
    int x_max = 0, y_max = 0, z_max = 0;
    int x_min = 0, y_min = 0, z_min = 0;

public:
    App_Stats() {}

    template<typename T>
    void gather_stats(QMatrix1D<T> vec) {
        T min = INT_MAX, max = INT_MIN;
        double sum = 0;
        for(int i = 0; i < vec.size(); i++) {
            if(min >= vec[i]) {
                x_min = i;
                min = vec[i];
            }

            if(max <= vec[i]) {
                x_max = i;
                max = vec[i];
            }

            sum += vec[i];
        }

        median_val = sum / (vec.size() + 1);
        max_val = max; min_val = min;

        double disp = 0;
        for(int i = 0; i < vec.size(); i++) {
            disp += (vec[i] - median_val) * (vec[i] - median_val);
        }
        sko = disp;
        dispersion = sqrt(disp / (vec.size() + 1));
    }

    template<typename T>
    void gather_stats(QMatrix2D<T> mat) {
        T min = INT_MAX, max = INT_MIN;
        double sum = 0;
        for(int i = 0; i < mat.size(); i++) {
            for(int j = 0; j < mat[0].size(); j++) {

                if(min > mat[i][j]) {
                    x_min = i;
                    y_min = j;
                    min = mat[i][j];
                }

                if(max < mat[i][j]) {
                    x_max = i;
                    y_max = j;
                    max = mat[i][j];
                }

                sum += mat[i][j];
            }
        }

        median_val = sum / (mat.size() + 1) / (mat[0].size() + 1);
        max_val = max; min_val = min;

        double disp = 0;
        for(int i = 0; i < mat.size(); i++) {
            for(int j = 0; j < mat[0].size(); j++) {
                disp += (mat[i][j] - median_val) * (mat[i][j] - median_val);
            }
        }
        sko = disp;
        dispersion = sqrt(disp / (mat.size() + 1) / (mat[0].size() + 1));
    }

};

struct ProfileResult {
    QVector<QPointF> points;
    QVector<QPoint> indices;
};

struct ProfileParameters {
    int x0 = 0;
    int y0 = 0;
    double sigma = 1.0;
    double xGradient = 0.0;
    double yGradient = 0.0;
    int numSigma = 3;
};

enum class attribute : int
{
    isEdge = 1, isSelectedEdge = 2
};

/*
 * EN:
 *  finds minimum value of vector<vector<T>>
 *  for any type that has > operator
 * RU:
 *  nahodit minimal'noye znacheniye v vector<vector<T>>
 *  dlya lyubogo tipa, u kotorogo yest' operator >
 */
template <typename T>
T min(Matrix2D<T> mat){
    T res = mat[0][0];
    for(size_t i = 0; i < mat.size(); i++){
        for(size_t j = 0; j < mat.size(); j++){
            if(res > mat[i][j]) res = mat[i][j];
        }
    }
    return res;
}

/*
 * EN:
 *  finds minimum value between two T values
 * RU:
 *  nahodit minimal'noye znacheniye mezhdu dvumya znacheniyami tipa T
 */
template <typename T>
T min(T val1, T val2){
    return val1 < val2 ? val1 : val2;
}

/*
 * EN:
 *  finds maximum value of vector<vector<T>>
 *  for any type that has > operator
 * RU:
 *  nahodit maksimal'noye znacheniye v vector<vector<T>>
 *  dlya lyubogo tipa, u kotorogo yest' operator >
 */
template <typename T>
T max(Matrix2D<T> mat){
    T res = mat[0][0];
    for(size_t i = 0; i < mat.size(); i++){
        for(size_t j = 0; j < mat.size(); j++){
            if(mat[i][j] > res) res = mat[i][j];
        }
    }
    return res;
}

/*
 * EN:
 *  finds maximum value between two T values
 * RU:
 *  nahodit maksimal'noye znacheniye mezhdu dvumya znacheniyami tipa T
 */
template <typename T>
T max(T val1, T val2){
    return val1 > val2 ? val1 : val2;
}

template<typename T>
void clearMatrix2D(Matrix2D<T>& mat) {
    for(auto& row : mat) {
        row.clear();
    }
    mat.clear();
}


//*********************************************************************************************************
// Sctructures for edge refinement using PCR algorithm / Structuri dlya utochneniya granits po metodu PCR
//*********************************************************************************************************

struct RefinementParameters {
    // Image data / Dannye izobrazheniya
    Matrix2D<double> A;        // First blurred image / Pervoe razmytoe izobrazhenie
    Matrix2D<double> B01;       // Second blurred image / Vtoroe razmytoe izobrazhenie
    double ex, ey;
    int NX, NY;

    // Profile parameters / Parametry profilya
    double sigma0;              // First sigma / Pervaya sigma
    double sigma01;              // First sigma / Pervaya sigma
    double sigma1;              // First sigma / Pervaya sigma
    double sigma2;              // First sigma / Pervaya sigma
    int n_sigma;                // Number of sigma steps / Kolichestvo shagov sigma
    int n_myu;                  // Number of mu points / Kolichestvo tochek mu
    int NN;                     // New profile length / Novaya dlina profilya
    int otstup;                 // Search range / Diapozon poiska
};

struct RefinementResult {
    QPointF refinedPosition;    // Refined position / Utochnennaya poziciya
    double FFF1_1;              // Shift value / Velichina sdviga
    double FFF1_0;              // Scale value / Velichina masshtaba
    double residual;            // Final residual / Finalnaya nevyazka
    bool success;               // Whether refinement succeeded / Uspeshno li utochnenie
};

#endif // IMAGEPROCESSINGTYPES_H
