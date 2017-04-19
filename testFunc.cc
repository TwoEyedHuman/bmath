#include <iostream>
#include "bMath.h"

double f(double x);

int main() {
	std::cout << "root is " << newtonMethod(1.2, f) << std::endl;
	return 0;
}

double f(double x) {
	return 3*pow(x,2.0)-12;
}

