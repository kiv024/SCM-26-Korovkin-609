#include "eqsolver.hpp"
#include <iostream>
#include <cstdio>
#include <omp.h>

double phi(double x, double y){
    return sin(x + y);
}

int is_inside(double x, double y){
    if((0 < x && x < 1 && pow(y, 2) < x) || (std::abs(x-1)<1e-9 && -1 < y && y < 1)){
        return 1;
    }
    if(pow(y, 2) > x || x>1){
        return -1;
    }
    return 1;
}

int main(){
    double x_min = -0.1;
    double x_max = 1.1;
    double y_min = -1.1;
    double y_max = 1.1;
    size_t N_x = 400;
    size_t N_y = 600;
    size_t iteration_no = 0;
    double time_taken = 0;

    auto solver = EqSolver(phi, is_inside, x_min, x_max, y_min, y_max, N_x, N_y);

    double start = omp_get_wtime();
    auto result = solver.solve_straightforward(1e-6, iteration_no);
    // auto result = solver.gradient_solve(1e-6, iteration_no);
    double finish = omp_get_wtime();
    time_taken = finish - start; //in seconds

    const char* file_name = "./computation_result.txt";
    FILE* file = fopen(file_name, "w");
    std::fprintf(file, "%f %zu\n", time_taken, iteration_no);
    for(size_t i = 0; i < N_x + 1; i++){
        for(size_t j = 0; j < N_y; j++){
            std::fprintf(file, "%f ", result.get(i, j));
        }
        std::fprintf(file, "%f\n", result.get(i, N_y));
    }
    fclose(file);

    // std::cout << result;
    // std::cout << iteration_no << std::endl;
}