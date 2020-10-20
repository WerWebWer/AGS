#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "./func.h"

using namespace std;

double func(double x) {
	return x * x;
}

class Pos {
public:
	double x;
	double y;
	
	Pos() {
		x = 0;
		y = 0;
	}
	Pos(double _x, double _y) {
		x = _x;
		y = _y;
	}
	Pos(double _x) {
		x = _x;
		y = func(x);
	}
	Pos(const Pos& pos) {
		x = pos.x;
		y = pos.y;
	}
	~Pos() {};

	bool operator <(const Pos& pos) { return x < pos.x; }

	friend bool operator< (Pos const& l, Pos const& r) { return l.x < r.x; }
	friend bool operator> (Pos const& l, Pos const& r) { return l.x > r.x; }
	friend bool operator<= (Pos const& l, Pos const& r) { return l.x <= r.x; }
	friend bool operator>= (Pos const& l, Pos const& r) { return l.x >= r.x; }
	friend bool operator== (Pos const& l, Pos const& r) { return l.x == r.x; }
	friend bool operator!= (Pos const& l, Pos const& r) { return l.x != r.x; }
	friend std::ostream& operator<< (std::ostream& out, const Pos& pos) {
		out << '(' << pos.x << ", " << pos.y << ')';
		return out;
	}
};


int main(){
	double e = 1e-16;
	double r = 50;
	double MAX_cur = 1000;

	Pos cur;
	vector<Pos> a;
	a.push_back(-1);
	a.push_back(2);

	double m = 0;

	double cur_e = a[1].x - a[0].x;

	double R_max = 0;
	double R_pos = 1;
	int count = 0;

	while (cur_e > e && count < MAX_cur) {
		double M = abs((a[1].y - a[0].y) / (a[1].x - a[0].x));
		for (int i = 2; i < a.size(); i++) {
			double max;
			max = abs((a[i].y - a[i - 1].y) / (a[i].x - a[i - 1].x));
			if (M > max) M = max;
		}

		m = (M == 0) ? 1 : r * M;

		if (m == 0) m = 1;
		else m = r * m;

		R_max = m * (a[1].x - a[0].x) + (pow((a[1].y - a[0].y), 2) / (m * (a[1].x - a[0].x))) - 2 * (a[1].y + a[0].y);

		for (int i = 2; i < a.size(); i++) {
			double R = m * (a[i].x - a[i - 1].x) + (pow((a[i].y - a[i - 1].y), 2) / (m * (a[i].x - a[i - 1].x))) - 2 * (a[i].y + a[i - 1].y);

			if (R > R_max) {
				R_max = R;
				R_pos = i;
			}
		}

		cur_e = a[R_pos].x - a[R_pos - 1].x;

		cur = Pos((a[R_pos].x + a[R_pos - 1].x) / 2 - (a[R_pos].y - a[R_pos - 1].y) / (2 * m));

		a.push_back(cur);
		sort(a.begin(), a.end());

		int pos = find(a.begin(), a.end(), cur) - a.begin();
		std::cout << "pos " << pos << " = " << a[pos] << std::endl;
		count++;
	}

	std::cout << "min = " << cur << std::endl;
}