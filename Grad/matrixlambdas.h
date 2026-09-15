#ifndef MATRIXLAMBDAS_H
#define MATRIXLAMBDAS_H

#include <cmath>
#include <limits>
#include <type_traits>
#include <algorithm>

/*
 * EN: Collection of functors for element-wise matrix operations
 * RU: Kollektsiya funktorov dlya poelementnykh matrichnykh operatsiy
 */

namespace MatrixLambdas {

// ============================================================================
// Arithmetic Operations / Aritmeticheskiye operatsii
// ============================================================================

/**
 * EN: Subtraction operator (a - b)
 * RU: Operatsiya vychitaniya (a - b)
 */
template<typename T>
struct Subtract {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const { return a - b; }
};

/**
 * EN: Subtraction with zero clamping (negative values become 0)
 * RU: Vychitaniye s ogranicheniyem do nulya (otritsatel'nyye znacheniya stanovatsya 0)
 */
template<typename T>
struct SubtractWithZeroClamp {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const {
        T res = a - b;
        return (res > 0) ? res : static_cast<T>(0);
    }
};

/**
 * EN: Addition operator (a + b)
 * RU: Operatsiya slozheniya (a + b)
 */
template<typename T>
struct Add {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const { return a + b; }
};

/**
 * EN: Addition with clamping to 255 (for image processing)
 * RU: Slozheniye s ogranicheniyem do 255 (dlya obrabotki izobrazheniy)
 */
template<typename T>
struct AddWithMaxClamp {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const {
        T res = a + b;
        return (res < 255) ? res : static_cast<T>(255);
    }
};

/**
 * EN: Multiplication operator (a * b)
 * RU: Operatsiya umnozheniya (a * b)
 */
template<typename T>
struct Multiply {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const { return a * b; }
};

/**
 * EN: Division operator with protection against division by zero
 * RU: Operatsiya deleniya s zashchitoy ot deleniya na nol'
 */
template<typename T>
struct Divide {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const {
        const T epsilon = static_cast<T>(1e-10);
        if (std::abs(b) < epsilon) {
            return std::numeric_limits<T>::max();
        }
        return a / b;
    }
};

// ============================================================================
// Comparison Operations / Operatsii sravneniya
// ============================================================================

/**
 * EN: Minimum of two values
 * RU: Minimum iz dvukh znacheniy
 */
template<typename T>
struct Min {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const { return std::min(a, b); }
};

/**
 * EN: Maximum of two values
 * RU: Maximum iz dvukh znacheniy
 */
template<typename T>
struct Max {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const { return std::max(a, b); }
};

// ============================================================================
// Bitwise Operations (for integer types only) / Pobitovyye operatsii (tol'ko dlya tselochislennykh tipov)
// ============================================================================

/**
 * EN: Bitwise AND (converts to int for operation)
 * RU: Pobitovoye I (preobrazuyet v int dlya operatsii)
 */
template<typename T>
struct BitwiseAND {
    T operator()(T a, T b) const {
        return static_cast<T>(static_cast<int>(a) & static_cast<int>(b));
    }
};

/**
 * EN: Bitwise OR (converts to int for operation)
 * RU: Pobitovoye ILI (preobrazuyet v int dlya operatsii)
 */
template<typename T>
struct BitwiseOR {
    T operator()(T a, T b) const {
        return static_cast<T>(static_cast<int>(a) | static_cast<int>(b));
    }
};

/**
 * EN: Bitwise XOR (converts to int for operation)
 * RU: Pobitovoye isklyuchayushcheye ILI (preobrazuyet v int dlya operatsii)
 */
template<typename T>
struct BitwiseXOR {
    T operator()(T a, T b) const {
        return static_cast<T>(static_cast<int>(a) ^ static_cast<int>(b));
    }
};

// ============================================================================
// Statistical Operations / Statisticheskiye operatsii
// ============================================================================

/**
 * EN: Average of two values
 * RU: Sredneye znacheniye dvukh chisel
 */
template<typename T>
struct Average {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const {
        return (a + b) / static_cast<T>(2);
    }
};

/**
 * EN: Absolute difference between two values
 * RU: Modul' raznosti mezhdu dvumya znacheniyami
 */
template<typename T>
struct Difference {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic type");
    constexpr T operator()(T a, T b) const {
        return std::abs(a - b);
    }
};

// ============================================================================
// Special Operations for Image Processing / Spetsial'nyye operatsii dlya obrabotki izobrazheniy
// ============================================================================

/**
 * EN: Weighted sum (for blending) - a * alpha + b * (1 - alpha)
 * RU: Vzveshennaya summa (dlya smeshivaniya) - a * alpha + b * (1 - alpha)
 */
template<typename T>
struct WeightedSum {
    double alpha;

    explicit WeightedSum(double alpha_ = 0.5) : alpha(std::max(0.0, std::min(1.0, alpha_))) {}

    T operator()(T a, T b) const {
        return static_cast<T>(a * alpha + b * (1.0 - alpha));
    }
};

/**
 * EN: Square of difference (for error calculation)
 * RU: Kvadrat raznosti (dlya rascheta oshibki)
 */
template<typename T>
struct SquaredDifference {
    T operator()(T a, T b) const {
        T diff = a - b;
        return diff * diff;
    }
};

} // namespace MatrixLambdas

#endif // MATRIXLAMBDAS_H
