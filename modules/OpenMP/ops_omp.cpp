// Copyright 2021 Smirnov Aleksandr
#include <omp.h>
#include "stdio.h"
#include "math.h"
#include "time.h"

#include <iostream>
#include <algorithm>
#include <limits>

#include "vector"
#include "iterator"
#include "../../modules/OpenMP/ops_omp.h"

#define NUM_THREADS 4

double ParallelOperations(double inter[2], double fun(double x), double r, double e) {
    double start = omp_get_wtime();
	
    omp_set_num_threads(NUM_THREADS);
    double left = inter[0];
    double right = inter[1];
    double min = left;
	
	double h = (right - left)/(double)NUM_THREADS;
	
	#pragma omp parallel
    {
		// nthreads= omp_get_num_threads();
        double min_local;
        #pragma omp for nowait
        for (int i = 0; i < NUM_THREADS; i++) {
            int thread_num = omp_get_thread_num();
            double lr[2] = { left + h * thread_num, left + h * (thread_num + 1) };
            min_local = AGS(lr, fun, r, e);
        }
        #pragma omp critical 
        {
            int thread_num = omp_get_thread_num();
            // std::cout << "!!!   " << thread_num << " x= " << min_local << " y= " << fun(min_local) << "\n";
            if (fun(min_local) <= fun(min)) {
                min = min_local;
            }
        }
    }
    
	std::cout << min << std::endl;
    double finish = omp_get_wtime();
    std::cout << "Time in OpenMP: " << finish - start << std::endl;
    return min;
}


double Rmax(std::vector<double> a) {
    double R = a[0];
    int count = 0;
    for (int i = 1; i < a.size(); i++) {
        if (a[i] > R) {
            R = a[i];
            count = i;
        }
    }
    return count;
};

double AGS(double inter[2], double fun(double x), double r, double e) {
    double left = inter[0];
    double right = inter[1];
    // std::cout << "[ " << left << " ; " << right << " ]" << std::endl;
    // double r = 1; //
    double M; //
    // double e = 0.000001; //
    double new_point;
    double m;
    int k = 2;
    double r_max = 0;
    std::vector <std::pair<double, double>> point;
    std::vector <double> R;
    std::pair<double, double> p1(left, fun(left));
    std::pair<double, double> p2(right, fun(right));
    point.insert(point.begin(), p1);
    point.push_back(p2);
    while ((point[r_max + 1].first - point[r_max].first) > e) {
        M = 0;
        for (int i = 0; i < k - 1; i++) {
            int tmpM = abs((point[i + 1].second - point[i].second) / (point[i + 1].first - point[i].first));
            if (M < tmpM) M = tmpM;
        }
        m = (M == 0) ? 1 : r * M;
        for (int i = 0; i < k - 1; i++) {
            R.push_back(m * (point[i + 1].first - point[i].first) + ((point[i + 1].second - point[i].second) * (point[i + 1].second - point[i].second)) / (m * (point[i + 1].first - point[i].first)) - 2 * (point[i + 1].second + point[i].second));
        }
        r_max = Rmax(R);
        R.clear();
        new_point = (point[r_max + 1].first + point[r_max].first) / 2 - (point[r_max + 1].second - point[r_max].second) / (2 * m);
        point.push_back(std::pair<double, double>(new_point, fun(new_point)));
        sort(point.begin(), point.end());
        k++;
    }
    // std::cout << point[r_max].first << " " << point[r_max].second << std::endl;
    return point[r_max].first;
}

double ParallelAGS(double inter[2], double fun(double x), double r, double e) {
    double left = inter[0];
    double right = inter[1];
    // std::cout << "[ " << left << " ; " << right << " ]" << std::endl;
    // double r = 1; //
    double M; //
    // double e = 0.000001; //
    double new_point;
    double m;
    int k = 2;
    double r_max = 0;
    std::vector <std::pair<double, double>> point;
    std::vector <double> R;
    std::pair<double, double> p1(left, fun(left));
    std::pair<double, double> p2(right, fun(right));
    point.insert(point.begin(), p1);
    point.push_back(p2);
    while ((point[r_max + 1].first - point[r_max].first) > e) {
        M = 0;
        #pragma omp parallel for reduction(max:M)
        for (int i = 0; i < k - 1; i++) {
            int tmpM = abs((point[i + 1].second - point[i].second) / (point[i + 1].first - point[i].first));
            if (M < tmpM) M = tmpM;
        }
        m = (M == 0) ? 1 : r * M;
        for (int i = 0; i < k - 1; i++) {
            R.push_back(m * (point[i + 1].first - point[i].first) + ((point[i + 1].second - point[i].second) * (point[i + 1].second - point[i].second)) / (m * (point[i + 1].first - point[i].first)) - 2 * (point[i + 1].second + point[i].second));
        }
        r_max = Rmax(R);
        R.clear();
        new_point = (point[r_max + 1].first + point[r_max].first) / 2 - (point[r_max + 1].second - point[r_max].second) / (2 * m);
        point.push_back(std::pair<double, double>(new_point, fun(new_point)));
        sort(point.begin(), point.end());
        k++;
    }
    // std::cout << point[r_max].first << " " << point[r_max].second << std::endl;
    return point[r_max].first;
}