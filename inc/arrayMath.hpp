/**
 * @file arrayMath.hpp
 * @author your name (you@domain.com)
 * @brief mathematical functions for arrays (len 3)
 * @version 0.1
 * @date 2022-03-19
 * 
 * @copyright Don't care Copyright (c) 2022
 * 
 */

#ifndef ARRAY_MATH_HPP
#define ARRAY_MATH_HPP

#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <sstream>

/*
// subtract array b from array a and return result
template <typename T, size_t size>
inline T* add(const T (&arrayA)[size], const T (&arrayB)[size]) {
    T* sum = new T[size];
    for (size_t i = 0; i < size; i++) {
        sum[i] = arrayA[i] + arrayB[i];
    }
    return sum;
}*/
/*
// subtract array b from array a and return result
template <typename T, size_t size>
inline T* add(const T* arrayA, const T* arrayB) {
    T* sum = new T[size];
    for (size_t i = 0; i < size; i++) {
        sum[i] = arrayA[i] + arrayB[i];
    }
    return sum;
}

template <typename T>
inline T* add3(const T* arrayA, const T* arrayB) {return add<T, 3>(arrayA, arrayB);} */

/**
 * @brief add array b to array a inplace (array a) (vector addition)
 * 
 * @tparam T type
 * @return T(&)[3] array of T of size 3
 */
template <typename T>
inline T (&add3Ip(T (&arrayA)[3], const T (&arrayB)[3]))[3] {
    for (size_t i = 0; i < 3; i++) {
        arrayA[i] = arrayA[i] + arrayB[i];
    }
    return arrayA;
}

/**
 * @brief add array a and array b to result (vector addition)
 * 
 * @tparam T 
 * @param arrayA 
 * @param arrayB 
 * @param result 
 */
template <typename T>
inline void add3(const T (&arrayA)[3], const T (&arrayB)[3], const T (&result)[3]) {
    for (size_t i = 0; i < 3; i++) {
        result[i] = arrayA[i] + arrayB[i];
    }
}

// subtracts b from a in a (Ip: inplace) -> to make this memsafe -> can be entered in another fct without generating a mem leak
// when arrayA and arrayB are well defined
/**
 * @brief subtract array b from array a inplace (array a) (vector subtraction)
 * 
 * @tparam T 
 * @return T(&)[3] 
 */
template <typename T>
inline T (&sub3Ip(T (&arrayA)[3], const T (&arrayB)[3]))[3] {
    // check this one out -> normally one new needs one delete
    //auto sum = std::make_unique<T[]>(3);
    for (size_t i = 0; i < 3; i++) {
        arrayA[i] = arrayA[i] - arrayB[i];
    }
    return arrayA;
}

// subtracts b from a in a (Ip: inplace) -> to make this memsafe -> can be entered in another fct without generating a mem leak
// when arrayA and arrayB are well defined
/**
 * @brief subtract array b from array a inplace (array a) (vector subtraction)
 * 
 * @tparam T 
 * @param arrayA 
 * @param arrayB 
 * @return T* 
 */
template <typename T>
inline T* sub3Ip(T* arrayA, const T* arrayB) {
    // check this one out -> normally one new needs one delete
    //auto sum = std::make_unique<T[]>(3);
    for (size_t i = 0; i < 3; i++) {
        arrayA[i] = arrayA[i] - arrayB[i];
    }
    return arrayA;
}

/**
 * @brief multiply by scalar in place return = arrayA
 * 
 * @tparam T 
 * @return T(&)[3] 
 */
template <typename T>
inline T (&mul3Ip(T (&arrayA)[3], T value))[3] {
    // check this one out -> normally one new needs one delete
    //auto sum = std::make_unique<T[]>(3);
    for (size_t i = 0; i < 3; i++) {
        arrayA[i] = arrayA[i] * value;
    }
    return arrayA;
}

/**
 * @brief multiply by scalar to result array
 * 
 * @tparam T 
 * @param arrayA 
 * @param value 
 * @param result 
 */
template <typename T>
inline void mul3(T (&arrayA)[3], T value, T (&result)[3]) {
    // check this one out -> normally one new needs one delete
    //auto sum = std::make_unique<T[]>(3);
    for (size_t i = 0; i < 3; i++) {
        result[i] = arrayA[i] * value;
    }
}

/*
// subtract array b from array a and return result
template <typename T, size_t size>
inline T* sub(const T (&arrayA)[size], const T (&arrayB)[size]) {
    T* sum = new T[size];
    for (size_t i = 0; i < size; i++) {
        sum[i] = arrayA[i] - arrayB[i];
    }
    return sum;
}*/

//template <typename T>
//inline T* sub3(const T* arrA, const T* arrB) {return sub<T,3>(arrA, arrB);}

/**
 * @brief euclidean norm / magnitude -> sqrt(sum(xi²)) return:  >= 0
 * 
 * @tparam T 
 * @param arrayA 
 * @return T 
 */
template <typename T>
inline T eNorm3(const T* arrayA) {
    T sum = 0;
    for (size_t i = 0; i < 3; i++) {
        sum += arrayA[i] * arrayA[i];
    }
    return pow(sum, .5);
}

// calculate the difference for an array of bounds and assign it to the first 3 elements of the input array
// use with caution
/**
 * @brief calculate the difference for an array of bounds and assign it to the first 3 elements of the input array use with caution
 * 
 * @tparam T 
 * @return T(&)[6] 
 */
template <typename T>
inline T (&bounds2DiffIp(T (&arr)[6]))[6] {
    arr[0] = arr[1] - arr[0];
    arr[1] = arr[3] - arr[2];
    arr[2] = arr[5] - arr[4];
    return arr;
}

/**
 * @brief assumes array length of 6 uses euclidean norm (eNorm3)
 * 
 * @tparam T 
 * @param arr 
 * @return T 
 */
template <typename T>
inline T bounds2DiagLen(T* arr){
    double diffBnds[3] = {arr[1]-arr[0], arr[3]-arr[2], arr[5]-arr[4]};
    return eNorm3(diffBnds);
}

/**
 * @brief convert an array to a string with " :: " as separator
 * @example array2Str([1,2,3,4]) -> 1 :: 2 :: 3 :: 4
 * 
 * @tparam T 
 * @tparam size 
 * @param arr 
 * @return std::string 
 */
template <typename T, size_t size>
inline std::string array2Str(T* arr) {
    std::stringstream ss;
    for (size_t i = 0; i < size-1; i++) {
        ss << arr[i] << " :: ";
    }
    ss << arr[size-1];
    return ss.str();
}

inline std::string array2Str3d(double* arr) {
    return array2Str<double, 3>(arr);
}

inline std::string array2Str6d(double* arr) {
    return array2Str<double, 6>(arr);
}

#endif // ARRAY_MATH_HPP