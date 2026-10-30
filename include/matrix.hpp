#pragma once

#include <iostream>
#include <random>
#include <vector>
#include <cmath>
  
template <typename T>
class Matrix{
private:
    size_t rows;
    size_t columns;
    
public:

    std::vector<T> data;
    
    Matrix() = default;

    size_t get_rows() const noexcept{
        return this->rows;
    }
    size_t get_columns() const noexcept{
        return this->columns;
    }
    
    Matrix(size_t rows, size_t columns){
        this->rows = rows;
        this->columns = columns;
        this->data.resize(rows*columns);
    }

    Matrix(std::vector<T> data, size_t rows, size_t columns){
        this->rows = rows;
        this->columns = columns;
        this->data.resize(rows*columns);
        #pragma omp parallel for proc_bind(spread)
        for(size_t i=0; i<rows*columns;i++){
            this->data[i] = data[i];
        }
    }

    Matrix(T* data_pointer, int data_size, size_t rows, size_t columns){
        this->rows = rows;
        this->columns = columns;
        this->data.resize(rows*columns);
        if(data_size != rows*columns){
            std::cerr<<"data_size != rows*columns"<<std::endl;
            std::exit(1);
        }
        #pragma omp parallel for proc_bind(spread)
        for(size_t i=0; i<rows*columns;i++){
            this->data[i] = data_pointer[i];
        }
    }

    const T get(size_t i, size_t j) const {
        return data[i * this->columns + j];
    }

    void set(size_t i, size_t j, T value){
        data[i * this->columns + j] = value;
    }
    
    // const T& operator[] (size_t i, size_t j) const{
    //     return data[i * this->columns + j];
    // }

    // T& operator[] (size_t i, size_t j){
    //     return data[i * this->columns + j];
    // }

    // Matrix<T>& operator= (const Matrix<T>& other){
    //     this->data = other.data;
    //     return *this;
    // }
};

template <typename T>
Matrix<T> operator+ (const Matrix<T>& mat1, const Matrix<T>& mat2){
    if ((mat1.get_rows() != mat2.get_rows()) || (mat1.get_columns() != mat2.get_columns())){
        std::cout << "uncorrelated sizes for +" << std::endl;
        exit(1);
    }
    auto result = Matrix<T>(mat1.get_rows(), mat1.get_columns());
    #pragma omp parallel for proc_bind(spread)
    for (size_t i = 0; i < mat1.get_rows(); i++){
        for (size_t j = 0; j < mat1.get_columns(); j++){
            result.set(i, j, mat1.get(i, j) + mat2.get(i, j));
        }
    }
    return result;
}

template <typename T>
Matrix<T> operator- (const Matrix<T>& mat1, const Matrix<T>& mat2){
    if ((mat1.get_rows() != mat2.get_rows()) || (mat1.get_columns() != mat2.get_columns())){
        std::cout << "uncorrelated sizes for -" << std::endl;
        exit(1);
    }
    Matrix<T> result = Matrix<T>(mat1.get_rows(), mat1.get_columns());
    #pragma omp parallel for proc_bind(spread)
    for (size_t i = 0; i < mat1.get_rows(); i++){
        for (size_t j = 0; j < mat1.get_columns(); j++){
            result.set(i, j, mat1.get(i, j) - mat2.get(i, j));
        }
    }
    return result;
}

template <typename T>
Matrix<T> operator* (T& scalar, const Matrix<T>& mat){
    Matrix<T> result = Matrix<T>(mat.get_rows(), mat.get_columns());
    #pragma omp parallel for proc_bind(spread)
    for (size_t i = 0; i < mat.get_rows(); i++){
        for (size_t j = 0; j < mat.get_columns(); j++){
            result.set(i, j, scalar * mat.get(i, j));
        }
    }
    return result;
}

template <typename T>
Matrix<T> operator* (const Matrix<T>& mat, T& scalar){
    Matrix<T> result = scalar*mat;
    return result;
}

template <typename T>
Matrix<T> operator* (const Matrix<T>& mat1, const Matrix<T>& mat2){
    if ((mat1.get_columns()) != (mat2.get_rows())){
        std::cout << "uncorrelated sizes for matrix multiplication" << std::endl;
        exit(1);
    }
    T temp;
    Matrix<T> result = Matrix<T>(mat1.get_rows(), mat2.get_columns());
    #pragma omp parallel for proc_bind(spread)
    for (size_t i = 0; i < mat1.get_rows(); i++){
        for (size_t j = 0; j < mat2.get_columns(); j++){
            temp = 0.0;
            for (size_t p = 0; p < mat1.get_columns(); p++){
                temp += mat1.get(i, p) * mat2.get(p, j);  
            }
            result.set(i, j, temp);
        }
    }
    return result;
}

// template<typename T>
// Matrix<T> fill_random(const Matrix<T>& mat, int random_lowest, int random_highest){
//     if (random_highest<random_lowest){
//         std::cout << "upper limit<lower limit" << std::endl;
//         exit(1);
//     }
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_int_distribution<int> dis(random_lowest, random_highest);
//     for (size_t i = 0; i < (mat.get_rows()); i++) {
//         for (size_t j = 0; j < (mat.get_columns()); j++){
//             mat[i, j] = dis(gen);
//         }
//     }
//     return mat;
// }

// template<typename T>
// std::ostream& operator<< (std::ostream& os, const Matrix<T>& mat){
//     for (size_t i = 0; i < mat.get_rows(); i++){
//         for(size_t j = 0; j < mat.get_columns(); j++){
//             os << mat[i, j] << ' ';
//         }
//         os << std::endl;
//     }
//     return os;
// }

