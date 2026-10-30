OMP :
	xlC_r -O3 -qstrict -qarch=pwr8 -qtune=pwr8 -qsmp=omp -I ./include -std=c++11 ./main.cpp ./source/eqsolver.cpp -o EqSolver

SEQ :
	xlC_r -O3 -qstrict -qarch=pwr8 -qtune=pwr8 -I ./include -std=c++11 ./main.cpp ./source/eqsolver.cpp -o EqSolver_seq