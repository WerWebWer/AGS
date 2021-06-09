// Copyright 2021 Smirnov Aleksandr
#include "stdio.h"
#include "math.h"
#include "time.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>

#include "../../modules/Sequence/func.h"

size_t Rmax(std::vector<double> *a){
    double *R = a->data();
    size_t count = 0;
    for (size_t i = 1; i < a->size(); i++) {
        if (*(a->data() + i) > * R) {
            R = a->data() + i;
            count = i;
        }
    }
    return count;
};

double GSA(double inter[2], double fun(double x), double r, double e){
    double left = inter[0];
    double right = inter[1];
    double M = 0;
    double new_point;
    double m; 
    size_t k = 2;
    size_t r_max = 0;
    std::vector <std::pair<double, double>> point{ std::pair<double, double>(left, fun(left)), 
                                                   std::pair<double, double>(right, fun(right)) };
    std::vector <double> R;
    while ((point[r_max + 1].first - point[r_max].first) > e) {
        M = 0;
        for (auto it1 = point.begin(), it2 = ++point.begin(); it2 != point.end(); it1++, it2++)
            M = fmax(M, abs((it2->second - it1->second) / (it2->first - it1->first)));
        m = (M == 0) ? 1 : r * M;
        for (auto it1 = point.begin(), it2 = ++point.begin(); it2 != point.end(); it1++, it2++)
            R.push_back(m * (it2->first - it1->first) + ((it2->second - it1->second) * (it2->second - it1->second)) / 
                        (m * (it2->first - it1->first)) - 2 * (it2->second + it1->second));
        r_max = Rmax(&R);
        R.clear();
        new_point = (point[r_max + 1].first + point[r_max].first) / 2 - (point[r_max + 1].second - point[r_max].second) / (2 * m);
        point.push_back(std::pair<double, double>(new_point, fun(new_point)));
        sort(point.begin(), point.end());
        k++;
    }
    return point[r_max].first;
}