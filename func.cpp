#include <cmath>
#include <algorithm>
#include <vector>
#include "./Pos.h"
#include "./func.h"

double function(double x) {
	return x * x;
}

double M_m(std::vector<Pos>& a, double r) {
	double M = abs((a[1].y - a[0].y) / (a[1].x - a[0].x));
	for (int i = 2; i < a.size(); i++) {
		double max;
		max = abs((a[i].y - a[i - 1].y) / (a[i].x - a[i - 1].x));
		if (M > max) M = max;
	}
	return (M == 0) ? 1 : r * M;
}

std::pair<int, double> func_R(std::vector<Pos>& a, double m) {
	int R_pos = 1;
	double R_max = m * (a[1].x - a[0].x) + (pow((a[1].y - a[0].y), 2) / (m * (a[1].x - a[0].x))) - 2 * (a[1].y + a[0].y);

	for (int i = 2; i < a.size(); i++) {
		double R = m * (a[i].x - a[i - 1].x) + (pow((a[i].y - a[i - 1].y), 2) / (m * (a[i].x - a[i - 1].x))) - 2 * (a[i].y + a[i - 1].y);

		if (R > R_max) {
			R_max = R;
			R_pos = i;
		}
	}
	return std::pair<int, double>(R_pos, R_max);
}