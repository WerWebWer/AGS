// Copyright 2021 Smirnov Aleksandr
#include <vector>
#include <string>
#include <utility>
#include <random>
#include <iostream>
#include "../../modules/Thread/ops_std.h"
#include "../../3rdparty/unapproved/unapproved.h"

std::mutex my_mutex;

size_t Rmax(std::vector<double>* a) {
    double* R = a->data();
    size_t count = 0;
    for (int i = 1; i < a->size(); i++) {
        if (*(a->data() + i) > * R) {
            R = a->data() + i;
            count = i;
        }
    }
    return count;
};

// true algorithm
void GSA(double inter[2], double (*fun)(double x), double r, double e, std::promise<int>&& pr) {
    double left = inter[0];
    double right = inter[1];
    double M = 0;
    double new_point;
    double m;
    size_t k = 2;
    size_t r_max = 0;
    std::vector <std::pair<double, double>> point{ std::pair<double, double>(left, fun(left)), std::pair<double, double>(right, fun(right)) };
    std::vector <double> R;
    while ((point[r_max + 1].first - point[r_max].first) > e) {
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
        sort(point.begin(), point.end());
        k++;
    }
    pr.set_value(round(point[r_max].first * 10));
};

// parallelization into segments
double ParallelGSA(double inter[2], double (*fun)(double x), double r, double e) {
    const int nthreads = std::thread::hardware_concurrency();
    const double delta = (inter[1] - inter[0]) / nthreads;
    std::promise<int>* promises = new std::promise<int>[nthreads];
    std::future<int>* futures = new std::future<int>[nthreads];
    std::thread* threads = new std::thread[nthreads];

    for (int i = 0; i < nthreads; i++) {
        futures[i] = promises[i].get_future();
        double lr[2] = { inter[0] + (double)i * delta,inter[0] + (double)(i + 1) * delta };
        threads[i] = std::thread(GSA, lr, fun, r, e, std::move(promises[i]));
        threads[i].join();
    }
    double min = inter[0];
    for (int i = 0; i < nthreads; i++) {
        double tmp = (double)futures[i].get() / 10.0;
        double fmin = fun(min);
        double ftmp = fun(tmp);
        if (ftmp <= fmin) {
            min = tmp;
        }
    }
    delete[]promises;
    delete[]futures;
    delete[]threads;
    return min;
};

void findM(std::vector<std::pair<double, double>>* point, size_t rank, std::promise<int>&& pr) {
    const int nthreads = std::thread::hardware_concurrency();
    std::pair<double, double>* point1 = point->data() + rank;
    std::pair<double, double>* point2 = point->data() + rank + 1;
    double M = 0;
    for (size_t i = rank + 1; i < point->size(); i += nthreads, point1 += nthreads, point2 += nthreads) {
        double tmpM = abs((point2->second - point1->second) / (point2->first - point1->first));
        if (M < tmpM) M = tmpM;
    }
    pr.set_value(round(M * 10000));
};

void calculateR(std::vector <double>* R, std::vector<std::pair<double, double>>* point, size_t rank, double m, std::promise<int>&& pr) {
    const int nthreads = std::thread::hardware_concurrency();
    std::pair<double, double>* point1 = point->data() + rank;
    std::pair<double, double>* point2 = point->data() + rank + 1;
    double* Ri = R->data() + rank;
    for (size_t i = rank; i < R->size(); i += nthreads, Ri += nthreads, point1 += nthreads, point2 += nthreads) {
        double tmp = (m * (point2->first - point1->first) + ((point2->second - point1->second) * (point2->second - point1->second)) / (m * (point2->first - point1->first)) - 2 * (point2->second + point1->second));
        *Ri = tmp;
    }
    pr.set_value(1);
};

// parallelization internals
double ParallelOperations(double inter[2], double (*fun)(double x), double r, double e) {
    const int nthreads = std::thread::hardware_concurrency();
    std::promise<int>* promisesM;
    std::future<int>* futuresM = new std::future<int>[nthreads];
    std::thread* threadsM = new std::thread[nthreads];

    std::promise<int>* promisesR;
    std::future<int>* futuresR = new std::future<int>[nthreads];
    std::thread* threadsR = new std::thread[nthreads];

    double left = inter[0];
    double right = inter[1];
    double M;
    double new_point;
    double m;
    size_t k = 2;
    int r_max = 0;
    std::vector <std::pair<double, double>> point{ std::pair<double, double>(left, fun(left)), std::pair<double, double>(right, fun(right)) };
    std::vector <double> R;
    while ((point[r_max + 1].first - point[r_max].first) > e) {
        promisesM = new std::promise<int>[nthreads];
        promisesR = new std::promise<int>[nthreads];
        M = 0;
        if (k >= nthreads + 1) {
            for (size_t i = 0; i < nthreads; i++) {
                futuresM[i] = promisesM[i].get_future();
                threadsM[i] = std::thread(findM, &point, i, std::move(promisesM[i]));
                threadsM[i].join();
            }
            for (size_t i = 0; i < nthreads && i < k - 1; i++) {
                int tmp = (double)futuresM[i].get() / 10000.0;
                if (M < tmp) M = tmp;
                futuresM[i].valid();
            }
        }
        else {
            for (size_t i = 0; i < k - 1; i++) {
                int tmpM = abs((point[i + 1].second - point[i].second) / (point[i + 1].first - point[i].first));
                if (M < tmpM) M = tmpM;
            }
        }

        m = (M == 0) ? 1 : r * M;

        R.push_back(0);

        if (k >= nthreads + 1) {
            int count = 0;
            for (size_t i = 0; i < nthreads; i++) {
                futuresR[i] = promisesR[i].get_future();
                threadsR[i] = std::thread(calculateR, &R, &point, i, m, std::move(promisesR[i]));
                threadsR[i].join();
                count += futuresR[i].get();
            }
        }
        else {
            for (size_t i = 0; i < k - 1; i++) {
                R[i] = m * (point[i + 1].first - point[i].first) + ((point[i + 1].second - point[i].second) * (point[i + 1].second - point[i].second)) / (m * (point[i + 1].first - point[i].first)) - 2 * (point[i + 1].second + point[i].second);
            }
        }

        r_max = Rmax(&R);
        new_point = (point[r_max + 1].first + point[r_max].first) / 2 - (point[r_max + 1].second - point[r_max].second) / (2 * m);
        size_t left = 0, right = point.size();
        while (left < right) {
            size_t mid = (left + right) / 2;
            if (new_point < point[mid].first)
                right = mid;
            else
                left = mid + 1;
        }
        point.insert(point.begin() + left, std::pair<double, double>(new_point, fun(new_point)));
        k++;
        delete[]promisesM;
        delete[]promisesR;
    }
    delete[]threadsM;
    delete[]futuresM;
    delete[]threadsR;
    delete[]futuresR;
    return point[r_max].first;
};

size_t* RmaxN(std::vector <std::pair<double, int>> R, size_t n) {
    size_t* r_max = new size_t[n];
    std::sort(R.begin(), R.end());
    for (size_t i = 0; i < n; i++)
        r_max[i] = R[R.size() - 1 - i].second;
    return r_max;
};

void calculateNewPoint(std::vector<std::pair<double, double>>* point, double (*fun)(double x), double *m, double *e, size_t *k, size_t *r, int i /*rank*/, double *fact, std::promise<int>&& pr) {
    double new_point = ((*(point->data() + *r + 1)).first + (*(point->data() + *r)).first) / 2 - ((*(point->data() + *r + 1)).second - (*(point->data() + *r)).second) / (2 * (*m));
    *(point->data() + *k + 1 + i) = std::pair<double, double>(new_point, fun(new_point));
    if ((*(point->data() + *r + 1)).first - (*(point->data() + *r)).first < *e) {
        pr.set_value(1);
        *fact = new_point;
    }
    else
        pr.set_value(0);
};

// parallel search for newpoints
double ParallelNewPoints(double inter[2], double (*fun)(double x), double r, double e) {
    const int nthreads = std::thread::hardware_concurrency();
    std::promise<int>* promises;
    std::future<int>* futures = new std::future<int>[nthreads];
    std::thread* threads = new std::thread[nthreads];
    double left = inter[0];
    double right = inter[1];
    double M;
    double new_point;
    double m;
    size_t k = 1;
    size_t* r_max;
    double fact = 0.;
    std::vector <std::pair<double, double>> point{ std::pair<double, double>(left, fun(left)), std::pair<double, double>(right, fun(right)) };
    std::vector <std::pair<double, int>> R = { std::pair<double, int>(0, 0) };
    bool flag = false;
    while (!flag) {
        promises = new std::promise<int>[nthreads];
        M = 0;
        for (auto it1 = point.begin(), it2 = ++point.begin(); it2 != point.end(); it1++, it2++)
            M = fmax(M, abs((it2->second - it1->second) / (it2->first - it1->first)));
        m = (M == 0) ? 1 : r * M;

        for (int i = 0; i < k; i++)
            R[i].first = m * (point[i + 1].first - point[i].first) + ((point[i + 1].second - point[i].second) * (point[i + 1].second - point[i].second)) / (m * (point[i + 1].first - point[i].first)) - 2 * (point[i + 1].second + point[i].second);

        size_t n = R.size() < nthreads ? R.size() : nthreads;
        for (size_t i = 0; i < n; i++)
            point.push_back(std::pair<double, double>(0, 0));
        r_max = RmaxN(R, n);
        for (int i = 0; i < n; i++) {
            futures[i] = promises[i].get_future();
            threads[i] = std::thread(calculateNewPoint, &point, fun, &m, &e, &k, &r_max[i], i, &fact, std::move(promises[i]));
            threads[i].join();
        }
        for (int i = 0; i < n; i++) {
            int tmp = futures[i].get();
            if (tmp == 1) {
                flag = true;
            }
            futures[i].valid();
        }
        sort(point.begin(), point.end());
        for (size_t i = 0; i < n; i++)
            R.push_back(std::pair<double, int>(0, R.size()));
        k += n;
        delete[]r_max;
        delete[]promises;
    }
    delete[]threads;
    delete[]futures;
    return fact;
}
