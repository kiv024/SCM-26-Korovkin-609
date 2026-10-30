#pragma once

#include <matrix.hpp>
#include <cmath>


class EqSolver{
public:
    EqSolver() = default;
    EqSolver(double (*phi)(double, double),int (*is_inside)(double, double), double x_min, double x_max, double y_min, double y_max, size_t N_x, size_t N_y);
    
    Matrix<double> solve_straightforward(double precision, size_t &iteration_no);
    Matrix<double> gradient_solve(double precision, size_t &iteration_no);
private:
    double h_x;
    double h_y;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    size_t N_x;
    size_t N_y;
    Matrix<double> Q;
    Matrix<double> S;
    Matrix<double> B;
    Matrix<double> D;
    double (*phi)(double, double);
    int (*is_inside)(double, double);

    void init_matricies();
    double sc_prod(Matrix<double>, Matrix<double>);
    Matrix<double> hadamard_inverse(Matrix<double>, Matrix<double>);
    Matrix<double> A(Matrix<double>);
};