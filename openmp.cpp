#include <omp.h>
#include <iostream>

#define NUM_THREADS 4

using namespace std;

int main(int argc, char* argv[]){
	
    omp_set_num_threads(NUM_THREADS);
	double val;
	double left = 1., right = 2.;
	double h = (right - left)/(double)NUM_THREADS;
	
#pragma omp parallel for reduction(min:val)
    {
		// nthreads= omp_get_num_threads();
		int thread_num = omp_get_thread_num();
        val = AGS(1+h*thread_num, 1+h*(thread_num + 1));
    }
    
	std::cout << val;

    return 0;
}