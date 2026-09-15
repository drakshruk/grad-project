#ifndef IMAGEPROCESSINGTYPES_H
#define IMAGEPROCESSINGTYPES_H

#include <QImage>
#include <QVector>
#include <QPoint>
#include <cmath>
#include <QDebug>
#include <vector>
#include <algorithm>
#include <limits>
#include <QFile>
#include <vector>
#include <algorithm>

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
    double mean_val = 0;      // Added: separate mean from median
    double dispersion = 0;
    double sko = 0;
    int x_max = 0, y_max = 0, z_max = 0;
    int x_min = 0, y_min = 0, z_min = 0;

public:
    // Default constructor is fine, but we can remove the explicit one
    // App_Stats() = default;  // Compiler-generated is sufficient

    template<typename T>
    void gather_stats(const QMatrix1D<T>& vec) {
        if (vec.empty()) return;

        T min = std::numeric_limits<T>::max();
        T max = std::numeric_limits<T>::min();
        double sum = 0.0;

        for(int i = 0; i < vec.size(); i++) {
            if(min >= vec[i]) {
                x_min = i;
                min = vec[i];
            }
            if(max <= vec[i]) {
                x_max = i;
                max = vec[i];
            }
            sum += static_cast<double>(vec[i]);
        }

        // Calculate mean correctly
        mean_val = sum / vec.size();
        max_val = static_cast<double>(max);
        min_val = static_cast<double>(min);

        // Calculate median - requires sorting a copy
        QVector<T> sorted = vec;
        std::sort(sorted.begin(), sorted.end());
        if (sorted.size() % 2 == 0) {
            median_val = static_cast<double>((sorted[sorted.size()/2 - 1] + sorted[sorted.size()/2]) / 2.0);
        } else {
            median_val = static_cast<double>(sorted[sorted.size()/2]);
        }

        // Calculate dispersion (variance) using mean, not median
        double disp = 0.0;
        for(int i = 0; i < vec.size(); i++) {
            double diff = static_cast<double>(vec[i]) - mean_val;
            disp += diff * diff;
        }
        dispersion = sqrt(disp / vec.size());
        sko = disp;  // Keeping as sum of squared differences
    }

    template<typename T>
    void gather_stats(const QMatrix2D<T>& mat) {
        if (mat.empty() || mat[0].empty()) return;

        T min = std::numeric_limits<T>::max();
        T max = std::numeric_limits<T>::min();
        double sum = 0.0;
        size_t rows = mat.size();
        size_t cols = mat[0].size();
        size_t totalPixels = rows * cols;

        for(size_t i = 0; i < rows; i++) {
            for(size_t j = 0; j < cols; j++) {
                T val = mat[i][j];
                if(min > val) {
                    x_min = static_cast<int>(i);
                    y_min = static_cast<int>(j);
                    min = val;
                }
                if(max < val) {
                    x_max = static_cast<int>(i);
                    y_max = static_cast<int>(j);
                    max = val;
                }
                sum += static_cast<double>(val);
            }
        }

        // Calculate mean correctly
        mean_val = sum / totalPixels;
        max_val = static_cast<double>(max);
        min_val = static_cast<double>(min);

        // Calculate median - flatten matrix into a vector
        QVector<T> flattened;
        flattened.reserve(totalPixels);
        for(size_t i = 0; i < rows; i++) {
            for(size_t j = 0; j < cols; j++) {
                flattened.push_back(mat[i][j]);
            }
        }
        std::sort(flattened.begin(), flattened.end());
        if (flattened.size() % 2 == 0) {
            median_val = static_cast<double>((flattened[flattened.size()/2 - 1] + flattened[flattened.size()/2]) / 2.0);
        } else {
            median_val = static_cast<double>(flattened[flattened.size()/2]);
        }

        // Calculate dispersion (variance) using mean, not median
        double disp = 0.0;
        for(size_t i = 0; i < rows; i++) {
            for(size_t j = 0; j < cols; j++) {
                double diff = static_cast<double>(mat[i][j]) - mean_val;
                disp += diff * diff;
            }
        }
        dispersion = sqrt(disp / totalPixels);
        sko = disp;  // Keeping as sum of squared differences
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
    isEdge = 1,
    isSelectedEdge = 2
};

/*
 * EN:
 *  finds minimum value of vector<vector<T>>
 *  for any type that has < operator
 * RU:
 *  nahodit minimal'noye znacheniye v vector<vector<T>>
 *  dlya lyubogo tipa, u kotorogo yest' operator <
 */
template <typename T>
T min(const Matrix2D<T>& mat){
    if (mat.empty() || mat[0].empty()) {
        return T{};
    }
    T res = mat[0][0];
    size_t rows = mat.size();
    size_t cols = mat[0].size();

    for(size_t i = 0; i < rows; i++){
        for(size_t j = 0; j < cols; j++){  // Fixed: use cols, not rows
            if(mat[i][j] < res) res = mat[i][j];
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
T max(const Matrix2D<T>& mat){
    if (mat.empty() || mat[0].empty()) {
        return T{};
    }
    T res = mat[0][0];
    size_t rows = mat.size();
    size_t cols = mat[0].size();

    for(size_t i = 0; i < rows; i++){
        for(size_t j = 0; j < cols; j++){  // Fixed: use cols, not rows
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



    // НОВЫЕ ПОЛЯ: предварительно отфильтрованные LoG изображения
    Matrix2D<double> A_log;       // LoG(A) — вычисляется ОДИН раз в test_001
    Matrix2D<double> B01_log;     // LoG(B01)


    double ex = 0.0, ey = 0.0;  // Gradient direction components / Komponenty napravleniya gradienta
    int NX = 0, NY = 0;          // Image dimensions / Razmery izobrazheniya

    // Profile parameters / Parametry profilya
    double sigma0 = 0.0;         // First sigma / Pervaya sigma
    double sigma01 = 0.0;        // Combined sigma for second blur / Kombinirovannaya sigma dlya vtorogo razmytiya
    double sigma1 = 0.0;         // First sigma for edge detection / Pervaya sigma dlya detekcii granic
    double sigma2 = 0.0;         // Second sigma for edge detection / Vtoraya sigma dlya detekcii granic
    int n_sigma = 0;             // Number of sigma steps / Kolichestvo shagov sigma
    int n_myu = 0;               // Number of mu points / Kolichestvo tochek mu
    int NN = 0;                  // New profile length / Novaya dlina profilya
    int otstup = 0;              // Search range / Diapozon poiska

    // Default constructor
    RefinementParameters() = default;
};

struct RefinementResult {
    QPointF refinedPosition;    // Refined position / Utochnennaya poziciya
    double FFF1_1 = 0.0;        // Shift value / Velichina sdviga
    double FFF1_0 = 0.0;        // Scale value / Velichina masshtaba
    double residual = 0.0;      // Final residual / Finalnaya nevyazka
    bool success = false;       // Whether refinement succeeded / Uspeshno li utochnenie

    // Default constructor
    RefinementResult() = default;
};

#endif // IMAGEPROCESSINGTYPES_H
