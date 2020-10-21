#include <vector>
#include <cmath>
#include <algorithm>
#include "./func.h"
#include "./Pos.h"

using namespace std;

int main(){
	//from user
	double e = 1e-10;
	double r = 1000;
	double MAX_count = 100;
	double left = -1;
	double right = 2;

	Pos cur;
	std::vector<Pos> a;
	int count = 0; // count iteration

	a.push_back(Pos(left, function(left)));
	a.push_back(Pos(right, function(right)));

	double cur_e = a[1].x - a[0].x;

	while (cur_e > e && count < MAX_count) {
		double M =       abs((a[1].y - a[1 - 1].y) / (a[1].x - a[1 - 1].x));
		for (int i = 2; i < a.size(); i++) {
			double max = abs((a[i].y - a[i - 1].y) / (a[i].x - a[i - 1].x));
			if (M > max) M = max;
		}

		double m = (M == 0) ? 1 : r * M;

		int R_pos = 1;
		double R_max = m * (a[1].x - a[1 - 1].x) + (pow((a[1].y - a[1 - 1].y), 2) / (m * (a[1].x - a[1 - 1].x))) - 2 * (a[1].y + a[1 - 1].y);
		for (int i = 2; i < a.size(); i++) {
			double R = m * (a[i].x - a[i - 1].x) + (pow((a[i].y - a[i - 1].y), 2) / (m * (a[i].x - a[i - 1].x))) - 2 * (a[i].y + a[i - 1].y);

			if (R > R_max) {
				R_max = R;
				R_pos = i;
			}
		}

		cur.x = (a[R_pos].x + a[R_pos - 1].x) / 2 - (a[R_pos].y - a[R_pos - 1].y) / (2 * m);
		cur.y = function(cur.x);

		a.push_back(cur);
		sort(a.begin(), a.end());
		if (count == 10)
			count += 0;
		int pos = find(a.begin(), a.end(), cur) - a.begin();
		std::cout << count << " pos " << pos << " = " << a[pos] << std::endl;
		cur_e = abs(a[R_pos].x - a[R_pos - 1].x);
		count++;
	}

	std::cout << std::endl << "min = " << cur << std::endl;
	Pos min = a[0];
	for (int i = 1; i < a.size(); i++) if (a[i].y < min.y) min = a[i];
	std::cout << "cur_min = " << min << std::endl;
}