#ifndef POS
#define POS

#include <iostream>

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
	Pos(const Pos& pos) {
		x = pos.x;
		y = pos.y;
	}
	~Pos() {};

	bool operator< (const Pos& pos) { return x < pos.x; }
	bool operator> (const Pos& pos) { return x > pos.x; }
	bool operator<= (const Pos& pos) { return x <= pos.x; }
	bool operator>= (const Pos& pos) { return x >= pos.x; }
	bool operator== (const Pos& pos) { return x == pos.x; }
	bool operator!= (const Pos& pos) { return x != pos.x; }
	friend std::ostream& operator<< (std::ostream& out, const Pos& pos) {
		out << '(' << pos.x << ", " << pos.y << ')';
		return out;
	}
};

#endif  // POS