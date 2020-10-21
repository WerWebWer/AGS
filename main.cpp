#include "./func.h"
#include "./Pos.h"
#include <vector>

using namespace std;

int main(){
	double e = 1e-16;
	double r = 50;
	double MAX_count = 1000;

	Pos cur;
	std::vector<Pos> a;

	a.push_back(Pos(-1, function(-1)));
	a.push_back(Pos(2, function(2)));

	double m = 0;

	double cur_e = a[1].x - a[0].x;

	double R_max = 0;
	int R_pos = 1;

	int count = 0;

	while (cur_e > e && count < MAX_count) {

		m = M_m(a, r); 

		//auto R = func_R(a, m);

		//R_pos = R.first;
		//R_max = R.second;

		R_max = m * (a[1].x - a[0].x) + (pow((a[1].y - a[0].y), 2) / (m * (a[1].x - a[0].x))) - 2 * (a[1].y + a[0].y);

		for (int i = 2; i < a.size(); i++) {
			double R = m * (a[i].x - a[i - 1].x) + (pow((a[i].y - a[i - 1].y), 2) / (m * (a[i].x - a[i - 1].x))) - 2 * (a[i].y + a[i - 1].y);

			if (R > R_max) {
				R_max = R;
				R_pos = i;
			}
		}
		//std::cout << "R_pos " << R_pos << " R_max " << R_max << std::endl;
		cur_e = a[R_pos].x - a[R_pos - 1].x;

		cur.x = (a[R_pos].x + a[R_pos - 1].x) / 2 - (a[R_pos].y - a[R_pos - 1].y) / (2 * m);
		cur.y = function(cur.x);

		a.push_back(cur);
		sort(a.begin(), a.end());

		int pos = find(a.begin(), a.end(), cur) - a.begin();
		//std::cout << "pos " << pos << " = " << a[pos] << std::endl;
		count++;
	}

	std::cout << "min = " << cur << std::endl;
}