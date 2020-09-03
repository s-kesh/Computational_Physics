#include<stdio.h>
#include<math.h>

#define NMAX 100
#define TOL 1e-8
#define X1 -2
#define X2 2

double absolute(double x)	{
	if (x < 0)
		return -1 * x;
	else
		return x;
}

double function(double x)	{
	return exp(-1 * x) + cos(x);
}

double average(double a, double b)	{
	return (a + b) / 2;
}

int main()	{
	double x1, x2, x;
	double f1, f2, fx;

	double absolute(double);
	double function(double);
	double average(double, double);

	printf("This program will find root of given function using Bisection Method\n");
	printf("Please enter lower bound and upper bound.\n");
	scanf("%f %f", &x1, &x2);

	f1 = function(x1);
	printf("%g\n", x1);
	printf("%g\n", f1);
	f2 = function(x2);
	printf("%g\n", x2);
	printf("%g\n", f2);

/****
	if (f1 * f2 < 0)	{
		for (int i = 0; i < NMAX; ++i)	{
			x = average(x1, x2);
			fx = function(x);
			if (fabs (fx) < TOL)
				break;
			if (function(x) * function(x1) < 0)
				x2 = x;
			else
				x1 = x;
		}
		printf("%g is root\n", x);
	}
	else if (function(x1) == 0)
		printf("%g is root\n", x1);
	else if (function(x2) == 0)
		printf("%g is root\n", x2);
	else
		printf("Select different interval");
****/

return 0;
}
