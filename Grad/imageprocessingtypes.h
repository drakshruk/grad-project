#ifndef IMAGEPROCESSINGTYPES_H
#define IMAGEPROCESSINGTYPES_H

#include <QImage>
#include <QVector>
#include <QPoint>
#include <cmath>
#include <QDebug>
#include <matrixlambdas.h>

template<typename T>
using Matrix2D = std::vector<std::vector<T>>;

const double pi = 4. * atan(1);

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

//template<typename T>
//void freeMatrix(T** matrix, int rows) {
//    if (matrix == nullptr) return;

//    for (int i = 0; i < rows; i++) {
//        if (matrix[i] != nullptr) {
//            delete[] matrix[i];
//            matrix[i] = nullptr;
//        }
//    }
//    delete[] matrix;
//}

//template<typename T>
//void copyMatrix(T** src, T** dest, int rows, int cols)
//{
//    for (int i = 0; i < rows; i++) {
//        for (int j = 0; j < cols; j++) {
//            dest[i][j] = src[i][j];
//        }
//    }
//}
template<typename T>
void clearMatrix2D(Matrix2D<T>& mat) {
    for(auto& row : mat) {
        row.clear();
    }
    mat.clear();
}

#endif // IMAGEPROCESSINGTYPES_H
