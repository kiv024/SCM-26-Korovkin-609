#include "eqsolver.hpp"

EqSolver::EqSolver(double (*phi)(double, double), int (*is_inside)(double, double), double x_min, double x_max, double y_min, double y_max, size_t N_x, size_t N_y) :
x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), N_x(N_x), N_y(N_y), phi(phi), is_inside(is_inside)
{
    this->h_x = (x_max-x_min)/N_x;
    this->h_y = (y_max-y_min)/N_y;
    init_matricies();

}


void EqSolver::init_matricies()
{
    this->S = Matrix<double>(N_x+1, N_y+1);
    this->Q = Matrix<double>(N_x+1, N_y+1);
    this->B = Matrix<double>(N_x+1, N_y+1);
    this->D = Matrix<double>(N_x+1, N_y+1);

    double area = this->h_x*this->h_y;
    #pragma omp parallel for shared(S, area) collapse(2)
    for (size_t i = 1; i <= this->N_x - 1; i++)
    {
        for (size_t j = 1; j <= this->N_y - 1; j++)
        {
            int sum = 0;

            double x_left = x_min + i * h_x - 0.5 * h_x;
            double x_right = x_min + i * h_x + 0.5 * h_x;
            double y_down = y_min + j * h_y - 0.5 * h_y;
            double y_up = y_min + j * h_y + 0.5 * h_y;

            auto in_ld = this->is_inside(x_left, y_down);
            auto in_lu = this->is_inside(x_left, y_up);
            auto in_ru = this->is_inside(x_right, y_up);
            auto in_rd = this->is_inside(x_right, y_down);

            sum += 1.0 + in_ld;
            sum += 1.0 + in_rd;
            sum += 1.0 + in_lu;
            sum += 1.0 + in_ru;

            this->S.set(i, j, (area/8.0)*sum);
        }
    }

    // #pragma omp parallel for
    // for(size_t i = 1; i <= N_x-1; i++){
    //     for(size_t j = 1; j <= N_y-1; j++){
    //         double x_l = x_min + i * h_x - 0.5 * h_x;
    //         double x_r = x_min + i * h_x + 0.5 * h_x;
    //         double y_d = y_min + j * h_y - 0.5 * h_y;
    //         double y_u = y_min + j * h_y + 0.5 * h_y;
    //         double area = h_x * h_y;
    //         int t_l = is_inside(x_l, y_u);
    //         int t_r = is_inside(x_r, y_u);
    //         int b_l = is_inside(x_l, y_d);
    //         int b_r = is_inside(x_r, y_d);
    //         if(t_l == 1 && t_r == 1 && b_l == 1 && b_r == 1){
    //             // S[i, j] = area;
    //             S.set(i, j, area);
    //             continue;
    //         }
    //         if(t_l == -1 && t_r == -1 && b_l == -1 && b_r == -1){
    //             // S[i, j] = 0;
    //             S.set(i, j, 0);
    //             continue;
    //         }
    //         if (t_l == 1 && t_r == 1 && b_l == -1 && b_r == 1){ //case 1
    //             // S[i, j] = area - 0.5 * std::abs((sqrt(x_l) + y_d) * (x_l - y_d*y_d));
    //             S.set(i, j, area - 0.5 * std::abs((sqrt(x_l) + y_d) * (x_l - y_d*y_d)));
    //             continue;
    //         }
    //         if(t_l == -1 && t_r == 1 && b_l == 1 && b_r == 1){ //case 2
    //             // S[i, j] = area - 0.5 * std::abs((sqrt(x_l) - y_u) * (x_l - y_u*y_u));
    //             S.set(i, j, area - 0.5 * std::abs((sqrt(x_l) - y_u) * (x_l - y_u*y_u)));
    //             continue;

    //         }
    //         if(t_l == -1 && t_r == 1 && b_l == -1 && b_r == -1){ //case 3 prim
    //             // S[i, j] = 0.5 * std::abs((sqrt(x_r) + y_u) * (x_r - y_u*y_u));
    //             S.set(i, j, 0.5 * std::abs((sqrt(x_r) + y_u) * (x_r - y_u*y_u)));
    //             continue;
    //         }
    //         if(t_l == -1 && t_r == -1 && b_l == -1 && b_r == 1){ //case 4 prim
    //             // S[i, j] = 0.5 * std::abs((sqrt(x_r) - y_d) * (x_r - y_d*y_d));
    //             S.set(i, j, 0.5 * std::abs((sqrt(x_r) - y_d) * (x_r - y_d*y_d)));
    //             continue;
    //         }
    //         if(t_l == -1 && t_r == -1 && b_l == 1 && b_r == 1){ //case 5
    //             // S[i, j] = area - 0.5 * (std::abs(sqrt(x_l) - y_u) + std::abs(sqrt(x_r) - y_u)) * std::abs(x_r - x_l);
    //             S.set(i, j, area - 0.5 * (std::abs(sqrt(x_l) - y_u) + std::abs(sqrt(x_r) - y_u)) * std::abs(x_r - x_l));
    //             continue;
    //         }
    //         if(t_l == 1 && t_r == 1 && b_l == -1 && b_r == -1){ //case 5 prim
    //             // S[i, j] = area - 0.5 * (std::abs(sqrt(x_l) + y_d) + std::abs(sqrt(x_r) + y_d)) * std::abs(x_r - x_l);
    //             S.set(i, j, area - 0.5 * (std::abs(sqrt(x_l) + y_d) + std::abs(sqrt(x_r) + y_d)) * std::abs(x_r - x_l));
    //             continue;
    //         }
    //         if(t_l == 1 && t_r == -1 && b_l == 1 && b_r == -1){ //case 6
    //             // S[i, j] = (1 - x_l) * std::abs(y_u - y_d);
    //             S.set(i, j, (1 - x_l) * std::abs(y_u - y_d));
    //             continue;
    //         };
    //         if(t_l == -1 && t_r == 1 && b_l == -1 && b_r == 1){ //case 6 prim
    //             // S[i, j] = 0.5 * (std::abs(x_r - y_u*y_u) + std::abs(x_r - y_d*y_d)) * std::abs(y_u - y_d);
    //             S.set(i, j, 0.5 * (std::abs(x_r - y_u*y_u) + std::abs(x_r - y_d*y_d)) * std::abs(y_u - y_d));
    //             continue;
    //         }
    //         if(t_l == 1 && t_r == -1 && b_l == -1 && b_r == -1){
    //             // S[i, j] = area;
    //             S.set(i, j, area);
    //             continue;
    //         }
    //         if(t_l == -1 && t_r == -1 && b_l == 1 && b_r == -1){
    //             // S[i, j] = area;
    //             S.set(i, j, area);
    //             continue;
    //         }
            

    //     }
    // }

    // std::cout << S;
    #pragma omp parallel for
    for(size_t i = 1; i <= N_x-1; i++){
        for(size_t j = 1; j <= N_y-1; j++){
            // Q[i, j] = (h_x * h_y - S[i, j]) / (h_x * h_y * pow(std::max(h_x, h_y), 4));
            // F[i, j] = 1;
            // B[i, j] = F[i, j] + Q[i, j] * phi(x_min + i * h_x, y_min + j * h_y);
            // D[i, j] = 2.0 / pow(h_x, 2) + 2.0 / pow(h_y, 2) + Q[i, j];
            Q.set(i, j, (h_x * h_y - S.get(i, j)) / (h_x * h_y * pow(std::max(h_x, h_y), 4)));
            B.set(i, j, 1 + Q.get(i, j) * phi(x_min + i * h_x, y_min + j * h_y));
            D.set(i, j, 2.0 / pow(h_x, 2) + 2.0 / pow(h_y, 2) + Q.get(i, j));
        }
    }
}

double EqSolver::sc_prod(Matrix<double> first, Matrix<double> second)
{
    double result = 0.0;
    #pragma omp parallel for reduction(+ : result)
    for(size_t i = 1; i <= N_x-1; i++){
        for(size_t j = 1; j <= N_y-1; j++){
            // result += h_x * h_y * first[i, j] * second[i, j];
            result += h_x * h_y * first.get(i, j) * second.get(i, j);
        }
    }   
    return result;
}

Matrix<double> EqSolver::hadamard_inverse(Matrix<double> multiplicand, Matrix<double> product)
{
    Matrix<double> result(N_x+1, N_y+1);
    #pragma omp parallel for
    for(size_t i = 1; i <= N_x-1; i++){
        for(size_t j = 1; j <=N_y-1; j++){
            // result[i, j] = product[i, j] / multiplicand[i, j];
            result.set(i, j, product.get(i, j) / multiplicand.get(i, j));
        }
    }   
    return result;

}

Matrix<double> EqSolver::A(Matrix<double> w)
{
    Matrix<double> result(N_x+1, N_y+1);
    #pragma omp parallel for
    for(size_t i = 1; i <= N_x-1; i++){
        for(size_t j = 1; j <= N_y-1; j++){
            // result[i, j] = (-w[i-1, j] + 2 * w[i, j] - w[i+1, j]) / pow(h_x, 2) + (-w[i, j-1] + 2 * w[i, j] - w[i, j+1]) / pow(h_y, 2) + Q[i, j] * w[i, j];
            result.set(i, j, (-w.get(i-1, j) + 2 * w.get(i, j) - w.get(i+1, j)) / pow(h_x, 2) + (-w.get(i, j-1) + 2 * w.get(i, j) - w.get(i, j+1)) / pow(h_y, 2) + Q.get(i, j) * w.get(i, j));
        }
    }   
    return result;
}

Matrix<double> EqSolver::gradient_solve(double precision, size_t &iteration_no)
{
    Matrix<double> result(N_x+1, N_y+1);
    Matrix<double> r_actual(N_x+1, N_y+1);
    Matrix<double> r_next(N_x+1, N_y+1);
    Matrix<double> w_actual(N_x+1, N_y+1);
    Matrix<double> w_next(N_x+1, N_y+1);
    Matrix<double> z_actual(N_x+1, N_y+1);
    Matrix<double> z_next(N_x+1, N_y+1);
    Matrix<double> p_actual(N_x+1, N_y+1);
    Matrix<double> p_next(N_x+1, N_y+1);
    double alpha;
    double beta;

    r_actual = B - A(w_actual);
    z_actual = hadamard_inverse(D, r_actual);
    p_actual = z_actual;
    alpha = sc_prod(z_actual, r_actual) / sc_prod(A(p_actual), p_actual);
    w_next = w_actual + alpha * p_actual;
    while(sc_prod(w_actual - w_next, w_actual - w_next) > 1e-16){
        w_actual = w_next;
        r_next = r_actual - alpha * A(p_actual);
        z_next = hadamard_inverse(D, r_next);
        beta = sc_prod(z_next, r_next) / sc_prod(z_actual, r_actual);
        p_next = z_next + beta * p_actual;
        alpha = sc_prod(z_next, r_next) / sc_prod(A(p_next), p_next);
        w_next = w_actual + alpha * p_next;

        p_actual = p_next;
        r_actual = r_next;
        z_actual = z_next;
    
        iteration_no++;
    }
    return w_next;
}

Matrix<double> EqSolver::solve_straightforward(double precision, size_t &iteration_no)
{
    Matrix<double> r(N_x+1, N_y+1);
    Matrix<double> w(N_x+1, N_y+1);
    Matrix<double> p(N_x+1, N_y+1);
    Matrix<double> A_applied_to_p(N_x+1, N_y+1);
    double alpha;
    double beta;

    double zr_next;
    
    double pAp;
    double pp;
    double weight = h_x * h_y;

    size_t max_iteration_no = (N_x - 1) * (N_y - 1);

    r = B - A(w);
    p = hadamard_inverse(D, r);
    pp = sc_prod(p, p);
    double zr_actual = sc_prod(p, r);
    A_applied_to_p = A(p);
    alpha = zr_actual / sc_prod(A_applied_to_p, p);
    double inv_h_x_squared = 1.0 / pow(h_x, 2);
    double inv_h_y_squared = 1.0 / pow(h_y, 2);
    while(pow(alpha, 2) * pp > 1e-16){
        zr_next = 0;
        #pragma omp parallel for reduction(+ : zr_next)
        for(size_t i = 1; i <= N_x-1; i++){
            for(size_t j = 1; j <=N_y-1; j++){
                w.set(i, j, w.get(i, j) + alpha * p.get(i, j));
                r.set(i, j, r.get(i, j) - alpha * A_applied_to_p.get(i, j));
                zr_next += r.get(i, j) / D.get(i, j) * r.get(i, j); 
            }
        }
        zr_next *= weight;

        beta = zr_next / zr_actual;
        zr_actual = zr_next;

        pAp = 0;
        pp = 0;

        #pragma omp parallel for
        for(size_t i = 1; i <= N_x-1; i++){
            for(size_t j = 1; j <= N_y-1; j++){
                p.set(i, j, r.get(i, j) / D.get(i, j) + beta * p.get(i, j));
            }
        }

        #pragma omp parallel for reduction(+ : pAp, pp)
        for(size_t i = 1; i <= N_x-1; i++){
            for(size_t j = 1; j <= N_y-1; j++){
                A_applied_to_p.set(i, j, (-p.get(i-1, j) - p.get(i+1, j)) * inv_h_x_squared + (-p.get(i, j-1) - p.get(i, j+1)) * inv_h_y_squared + D.get(i, j) * p.get(i, j));
                pAp += A_applied_to_p.get(i, j) * p.get(i, j); 
                pp += p.get(i, j) * p.get(i, j);
            }
        }
        pAp *= weight;
        pp *= weight;
        alpha = zr_next / pAp;

    
        iteration_no++;
        if(iteration_no > max_iteration_no){
            break;
        }
    }
    return w;
}
