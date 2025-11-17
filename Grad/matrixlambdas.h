#ifndef MATRIXLAMBDAS_H
#define MATRIXLAMBDAS_H

#include <cmath>
#include <limits>
#include <type_traits>

namespace MatrixLambdas {

// Use structs with operator() instead of template lambdas
template<typename T>
struct Subtract {
    constexpr T operator()(T a, T b) const { return a - b; }
};

template<typename T>
struct SubtractWithZeroClamp {
    constexpr T operator()(T a, T b) const {
        T res = a - b;
        return res > 0 ? res : static_cast<T>(0);
    }
};

template<typename T>
struct Add {
    constexpr T operator()(T a, T b) const { return a + b; }
};

template<typename T>
struct AddWithMaxClamp {
    constexpr T operator()(T a, T b) const {
        T res = a + b;
        return res < 255 ? res : static_cast<T>(255);
    }
};

template<typename T>
struct Multiply {
    constexpr T operator()(T a, T b) const { return a * b; }
};

template<typename T>
struct Divide {
    constexpr T operator()(T a, T b) const {
        return (b != 0 && std::abs(b) > 1e-6) ? a / b : std::numeric_limits<T>::max();
    }
};

template<typename T>
struct Min {
    constexpr T operator()(T a, T b) const { return a < b ? a : b; }
};

template<typename T>
struct Max {
    constexpr T operator()(T a, T b) const { return a > b ? a : b; }
};

// For bitwise operations, we need to specialize or handle type conversion
template<typename T>
struct BitwiseAND {
    constexpr T operator()(T a, T b) const {
        return static_cast<T>(static_cast<int>(a) & static_cast<int>(b));
    }
};

template<typename T>
struct BitwiseOR {
    constexpr T operator()(T a, T b) const {
        return static_cast<T>(static_cast<int>(a) | static_cast<int>(b));
    }
};

template<typename T>
struct BitwiseXOR {
    constexpr T operator()(T a, T b) const {
        return static_cast<T>(static_cast<int>(a) ^ static_cast<int>(b));
    }
};

template<typename T>
struct Average {
    constexpr T operator()(T a, T b) const { return (a + b) / static_cast<T>(2); }
};

template<typename T>
struct Difference {
    constexpr T operator()(T a, T b) const { return std::abs(a - b); }
};

} // namespace MatrixLambdas

#endif
