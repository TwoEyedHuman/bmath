#include <map>
#include <stdexcept>

const int ZERO_INT = 0;
const int MAX_INT = ~ZERO_INT;

double abs(double a);
void swap(double &a, double &b);
double bisectionMethod( double a, double b, double exp, double (*f)(double, double), double secondParam );
double bisectionMethod( double a, double b, double exp, double (*f)(double) );
double sqrt(double a);
double secantSlope(double a, double (*f)(double));
double secantSlope(double a, double (*f)(double), double epsilon);
double rootFinder(double a, double b);
//void rootFinder(double a, double b, double c, double & arr[2]);
double summation(int lowerLimit, int higherLimit, double (*f)(double));
double newtonMethod(double initEst, double (*f)(double));

double newtonMethod(double initEst, double (*f)(double) ){
	double x_n = initEst, epsilon = 0.01;
	while ( abs(f(x_n)) > epsilon) {
		x_n = x_n - (( f(x_n) )/( secantSlope(x_n, f) ));
	}
	return x_n;
}

double summation(int lowerLimit, int upperLimit, double (*f)(double)){
	double runsum = 0;
	for (int i=lowerLimit; i <=upperLimit; i++) {
		runsum = runsum + f(i);
	}
	return runsum;
}


double rootFinder(double a, double b) {
	if ( a == 0 ) {
		throw std::invalid_argument( "invalid coefficient: a" );
	}
	return -1*(b/a);
}

/*void rootFinder(double a, double b, double c, double & arr[2]) {
	if ( a == 0 ) {
		throw std::invalid_argument( "invalid coefficient: a" );
	}
	if ( pow(b,2) < 4*a*c ) {
		throw std:: std::invalid_argument( "invalid coefficients: b^2 < 4ac" );
	}
	arr[0] = (-1*b + sqrt(pow(b,2) - 4*a*c))/(2*a);
	arr[1] = (-1*b - sqrt(pow(b,2) - 4*a*c))/(2*a);
}
*/
double secantSlope(double a, double (*f)(double)){
	double epsilon = 0.01;
	return (f(a+epsilon/2) - f(a-epsilon/2))/epsilon;
}

double secantSlope(double a, double (*f)(double), double epsilon) {
	return (f(a+epsilon/2) - f(a-epsilon/2))/epsilon;
}

double abs(double a) {
	if ( a < 0 ) {
		return -1 * a;
	}
	else {
		return a;
	}
}

void swap(double &a, double &b) {
	double tmp = b;
	b = a;
	a = tmp;
}

double bisectionMethod( double a, double b, double exp, double (*f)(double, double), double secondParam ) {
	double epsilon = 0.00001;
	if (a > b) {
		swap(a,b);
	}
	double diff = b-a;
	while (diff > epsilon) {
		if (f(diff/2 + a, secondParam) > exp) {
			b = diff/2 + a;
		}
		else {
			a = diff/2 + a;
		}
		diff = b-a;
	}
	return a + diff;
}

double bisectionMethod( double a, double b, double exp, double (*f)(double) ) {
	double epsilon = 0.00001;
	if (a > b) {
		swap(a,b);
	}
	double diff = b-a;
	while (diff > epsilon) {
		if (f(diff/2 + a) > exp) {
			b = diff/2 + a;
		}
		else {
			a = diff/2 + a;
		}
		diff = b-a;
	}
	return a + diff;
}

template <class NumericType>
//NumericType pow(NumericType base, NumericType exponent);
NumericType pow(NumericType base, NumericType exponent) {
	NumericType runprod = base;
	for (int i=0; i<exponent-1; i++) {
		runprod = runprod * base;
	}
	return runprod;
}

template <class NumericType>
NumericType pow(NumericType base, NumericType exponent, bool overflowIndicate) {
	NumericType runprod = base;
	for (int i=0; i<exponent-1; i++) {
		runprod = runprod * base;
	}
	return runprod;
}

double sqrt(double a) {
	return bisectionMethod(0, 10000000000, a, pow, 2);
}
