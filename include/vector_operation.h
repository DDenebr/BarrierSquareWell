#ifndef VECTOR_OPERATION_INCLUDED
#define VECTOR_OPERATION_INCLUDED

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>

using std::array;

//Sign of variable
template <typename T> inline constexpr
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
};

//Sign of array
template <typename T, std::size_t size> inline constexpr
array<int, size> sgn_array(const array<T, size>& a){
    array<int, size> sign{};
    std::transform(a.begin(), a.end(), sign.begin(), sgn<T>);
    return sign;
};

//Sum of array
template <typename T, std::size_t size> inline constexpr
array<T, size> operator+(const array<T, size>& a, const array<T, size>& b){
    array<T, size> sum{};
    std::transform(a.begin(), a.end(), b.begin(), sum.begin(), std::plus<T>());
    return sum;
};
template <typename T, std::size_t size> inline constexpr
array<T, size> operator+(const array<T, size>& a, const T& b){
    array<T, size> sum{};
    std::transform(a.begin(), a.end(), sum.begin(), std::bind(std::plus<T>(), std::placeholders::_1, b));
    return sum;
};
template <typename T, std::size_t size> inline constexpr
array<T, size> operator+(const T& a, const array<T, size>& b){
    return b + a;
};

//Opposite of array
template <typename T, std::size_t size> inline constexpr
array<T, size> operator-(const array<T, size>& a){
    array<T, size> opposite{};
    std::transform(a.begin(), a.end(), opposite.begin(), std::negate<T>());
    return opposite;
};

//Difference of array
template <typename T, std::size_t size> inline constexpr
array<T, size> operator-(const array<T, size>& a, const array<T, size>& b){
    array<T, size> diff{};
    std::transform(a.begin(), a.end(), b.begin(), diff.begin(), std::minus<T>());
    return diff;
};
template <typename T, std::size_t size> inline constexpr
array<T, size> operator-(const array<T, size>& a, const T& b){
    array<T, size> diff{};
    std::transform(a.begin(), a.end(), diff.begin(), std::bind(std::minus<T>(), std::placeholders::_1, b));
    return diff;
};

//Product of array
template <typename T, std::size_t size> inline constexpr
array<T, size> operator*(const array<T, size>& a, const array<T, size>& b){
    array<T, size> prod{};
    std::transform(a.begin(), a.end(), b.begin(), prod.begin(), std::multiplies<T>());
    return prod;
};
template <typename T, std::size_t size> inline constexpr
array<T, size> operator*(const array<T, size>& a, const T& b){
    array<T, size> prod{};
    std::transform(a.begin(), a.end(), prod.begin(), std::bind(std::multiplies<T>(), std::placeholders::_1, b));
    return prod;
};
template <typename T, std::size_t size> inline constexpr
array<T, size> operator*(const T& a, const array<T, size>& b){
    return b * a;
};

//Quotient of array
template <typename T, std::size_t size> inline constexpr
array<T, size> operator/(const array<T, size>& a, const array<T, size>& b){
    array<T, size> quot{};
    std::transform(a.begin(), a.end(), b.begin(), quot.begin(), std::divides<T>());
    return quot;
};
template <typename T, std::size_t size> inline constexpr
array<T, size> operator/(const array<T, size>& a, const T& b){
    array<T, size> quot{};
    std::transform(a.begin(), a.end(), quot.begin(), std::bind(std::divides<T>(), std::placeholders::_1, b));
    return quot;
};

//Dot product
template <typename T, std::size_t size> inline constexpr
T dot(const array<T, size>& a, const array<T, size>& b){
    array<T, size> prod(a * b);
    return std::accumulate(prod.begin(), prod.end(), T(0));
};

template <typename T> inline constexpr
array<T, 3> cross(const array<T, 3>& a, const array<T, 3>& b){
    return array<T, 3>{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
};

//Return Euclidean norm (norm2)
template <typename T, std::size_t size> inline constexpr
T norm2(const array<T, size>& val){
    return std::sqrt(std::accumulate(val.begin(), val.end(), T(0), 
        [](const T& sum_of_square, const T& i){ return sum_of_square + i * i;}));
};

#endif