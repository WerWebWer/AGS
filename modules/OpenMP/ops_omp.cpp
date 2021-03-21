// Copyright 2021 Smirnov Aleksandr
#include <omp.h>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include "../../modules/OpenMP/ops_omp.h"

#define NUM_THREADS 4

int getParallelOperations(double left, double right) {
 //   double start = omp_get_wtime();
	//
 //   omp_set_num_threads(NUM_THREADS);
	//double val;
	//
	//double h = (right - left)/(double)NUM_THREADS;
	//
	//#pragma omp parallel for reduction(min:val)
 //   {
	//	// nthreads= omp_get_num_threads();
	//	int thread_num = omp_get_thread_num();
 //       val = AGS(1 + h * thread_num, 1 + h * (thread_num + 1));
 //   }
 //   
	//std::cout << val << std::endl;
 //   double finish = omp_get_wtime();
 //   std::cout << "How measure time in OpenMP: " << finish - start << std::endl;
 //   return reduction_elem;
    return 1;
}
