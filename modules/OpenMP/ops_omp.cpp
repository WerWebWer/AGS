// Copyright 2021 Smirnov Aleksandr
#include <omp.h>
#include "stdio.h"
#include "math.h"
#include "time.h"

#include <iostream>
#include <algorithm>
#include <limits>
#include <vector>
#include <iterator>

#include "../../modules/OpenMP/ops_omp.h"

#define NUM_THREADS 12

size_t Rmax(std::vector<double>* a) {
    double* R = a->data();
    size_t count = 0;
    for (size_t i = 1; i < a->size(); i++) {
        if (*(a->data() + i) > * R) {
            R = a->data() + i;
            count = i;
        }
    }
    return count;
};

// true algorithm
double GSA(double inter[2], double (*fun)(double x), double r, double e) {
    double left = inter[0];
    double right = inter[1];
    double M = 0;
    double new_point;
    double m;
    size_t k = 2;
    size_t r_max = 0;
    std::vector <std::pair<double, double>> point{ std::pair<double, double>(left, fun(left)), std::pair<double, double>(right, fun(right)) };
    std::vector <double> R;
    while (abs(point[r_max + 1].first - point[r_max].first) > e) {
        M = 0;
        for (auto it1 = point.begin(), it2 = ++point.begin(); it2 != point.end(); it1++, it2++)
            M = fmax(M, abs((it2->second - it1->second) / (it2->first - it1->first)));
        m = (M == 0) ? 1 : r * M;
        for (auto it1 = point.begin(), it2 = ++point.begin(); it2 != point.end(); it1++, it2++)
            R.push_back(m * (it2->first - it1->first) + ((it2->second - it1->second) * (it2->second - it1->second)) / (m * (it2->first - it1->first)) - 2 * (it2->second + it1->second));
        r_max = Rmax(&R);
        R.clear();
        new_point = (point[r_max + 1].first + point[r_max].first) / 2 - (point[r_max + 1].second - point[r_max].second) / (2 * m);
        point.push_back(std::pair<double, double>(new_point, fun(new_point)));
        std::sort(point.begin(), point.end());
        k++;
    }
    return point[r_max].first;
};

// parallelization into segments
double ParallelGSA(double inter[2], double (*fun)(double x), double r, double e) {
    omp_set_num_threads(NUM_THREADS);
    double left = inter[0];
    double right = inter[1];
    double min = left;

    double h = (right - left) / (double)NUM_THREADS;

#pragma omp parallel
    {
        double min_local;
#pragma omp for nowait
        for (int i = 0; i < NUM_THREADS; i++) {
            int thread_num = omp_get_thread_num();
            double lr[2] = { left + h * thread_num, left + h * (thread_num + 1) };
            min_local = GSA(lr, fun, r, e);
        }
#pragma omp critical 
        {
            if (fun(min_local) <= fun(min)) {
                min = min_local;
            }
        }
    }

    return min;
};

// parallelization internals
double ParallelOperations(double inter[2], double (*fun)(double x), double r, double e) {
    double left = inter[0];
    double right = inter[1];

    double M;
    double new_point;
    double m;
    int k = 2;
    int r_max = 0;
    std::vector <std::pair<double, double>> point{ std::pair<double, double>(left, fun(left)), std::pair<double, double>(right, fun(right)) };
    std::vector <double> R;
    while ((point[r_max + 1].first - point[r_max].first) > e) {
        M = 0;
#pragma omp parallel
        {
#pragma omp for nowait
            for (int i = 0; i < k - 1; i++) {
                double tmpM = abs((point[i + 1].second - point[i].second) / (point[i + 1].first - point[i].first));
#pragma omp critical
                if (M < tmpM) M = tmpM;
            }
        }
        m = (M == 0) ? 1 : r * M;

        R.push_back(0);

#pragma omp for  //nowait
        for (int i = 0; i < k - 1; i++)
            R[i] = m * (point[i + 1].first - point[i].first) + ((point[i + 1].second - point[i].second) * (point[i + 1].second - point[i].second)) / (m * (point[i + 1].first - point[i].first)) - 2 * (point[i + 1].second + point[i].second);

        r_max = Rmax(&R);
        new_point = (point[r_max + 1].first + point[r_max].first) / 2 - (point[r_max + 1].second - point[r_max].second) / (2 * m);
        point.push_back(std::pair<double, double>(new_point, fun(new_point)));
        std::sort(point.begin(), point.end());
        k++;
    }
    return point[r_max].first;
};

size_t* RmaxN(std::vector <std::pair<double, int>> R, size_t n) {
    size_t* r_max = new size_t[n];
    std::sort(R.begin(), R.end());
    for (size_t i = 0; i < n; i++)
        r_max[i] = R[R.size() - 1 - i].second;
    return r_max;
};

// parallel search for newpoints
double ParallelNewPoints(double inter[2], double (*fun)(double x), double r, double e) {
    double left = inter[0];
    double right = inter[1];
    double M;
    double new_point = 0;
    double m;
    int k = 1;
    size_t* r_max;
    std::vector <std::pair<double, double>> point{ std::pair<double, double>(left, fun(left)), std::pair<double, double>(right, fun(right)) };
    std::vector <std::pair<double, int>> R = { std::pair<double, int>(0, 0) };
    bool flag = false;
    while (!flag) {
        M = 0;
        for (auto it1 = point.begin(), it2 = ++point.begin(); it2 != point.end(); it1++, it2++)
            M = fmax(M, abs((it2->second - it1->second) / (it2->first - it1->first)));
        m = (M == 0) ? 1 : r * M;

        for (int i = 0; i < k; i++)
            R[i].first = m * (point[i + 1].first - point[i].first) + ((point[i + 1].second - point[i].second) * (point[i + 1].second - point[i].second)) / (m * (point[i + 1].first - point[i].first)) - 2 * (point[i + 1].second + point[i].second);

        int n = R.size() < NUM_THREADS ? R.size() : NUM_THREADS;
        omp_set_num_threads(n);
        for (int i = 0; i < n; i++)
            point.push_back(std::pair<double, double>(0, 0));
        r_max = RmaxN(R, n);
#pragma omp parallel
        {
            double new_point_local;
#pragma omp for
            for (int i = 0; i < n; i++) {
                size_t r = r_max[i];
                new_point_local = (point[r + 1].first + point[r].first) / 2 - (point[r + 1].second - point[r].second) / (2 * m);
                point[k + 1 + i] = std::pair<double, double>(new_point_local, fun(new_point_local));
#pragma omp critical
                if (point[r + 1].first - point[r].first < e) {
                    new_point = new_point_local;
                    flag = true;
                }
            }
        }
        std::sort(point.begin(), point.end());

        for (int i = 0; i < n; i++)
            R.push_back(std::pair<double, int>(0, R.size()));

        k += n;
        delete[]r_max;
    }
    return new_point;
};
