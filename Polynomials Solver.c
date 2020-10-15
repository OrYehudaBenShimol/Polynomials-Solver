#include<stdio.h>
#include<conio.h>
#include<math.h>

#define EPSILON 0.00001
#define PI 3.14159265

double BisectionMethod();				// find the real root with bisection method.
double SixthPolynom(double x);          //6th degree polynomial function.
double FifthPolynom(double x);          //5th degree polynomial function.
int solve_quadratic(double a, double b, double c, double *r1, double *i1, double *r2, double *i2);
void divide_cubic(double a, double *b, double *c, double x);
void cubic_root(double x, double y, double *r, double *i);
void square_root(double x, double y, double *r, double *i);
int solve_cubic(double a1, double a2, double a3, double a4, double *x1, double *y1, double *x2, double *y2, double *x3, double *y3);
void solve_biquadratic(double a, double b, double c, double *x1, double *y1, double *x2, double *y2, double *x3, double *y3, double *x4, double *y4);
int solve_quartic(double a4, double a3, double a2, double a1, double a0, double *x1, double *y1, double *x2, double *y2, double *x3, double *y3, double *x4, double *y4);
void complex_mult(double a, double b, double c, double d, double *x, double *y);
void complex_add(double a, double b,double c, double d, double *x, double *y);
void test_solution(double a, double b,double c, double d, double e, double x, double y,double *result1, double *result2);
double AX5, BX4, CX3, DX2, EX, F,GX6;
double A, B, C, D, E, F1, G;

int main()
{
	int i;
	double first_poly[6], poly_after_div[5], root1[2], reminder, derivative[7], minimum[5];
	double zero = 0, a=0, b=0, first_root = 0,temp1=0,temp2=0,temp3=0,temp4=0,temp5=0,min=0;
	double r, j;
	double a1, b1, c1, d1, e1;
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;
	int no_of_solutions;
	int flag;
	double result1, result2;

	printf("Please enter 6 Coefficients of the polynom (1st coefficient must be greater then zero):\n");
	scanf("%lf%lf%lf%lf%lf%lf", &AX5, &BX4, &CX3, &DX2, &EX, &F);

	while (AX5 <= 0) //if the first coefficient is equal to zero or less then zero , print the error massage below and ask for new input.
	{
		printf("Error, the first coefficient must be greater then zero.\n");
		printf("Please enter 6 Coefficients of the polynom (1st coefficient must be greater then zero):\n");
		scanf("%lf%lf%lf%lf%lf%lf", &AX5, &BX4, &CX3, &DX2, &EX, &F);
	}
	
	first_root=BisectionMethod();
	printf("\nBisection returned %.5lf\n", first_root);

	// put the coefficients into another array.
	first_poly[0] = AX5;
	first_poly[1] = BX4;
	first_poly[2] = CX3;
	first_poly[3] = DX2;
	first_poly[4] = EX;
	first_poly[5] = F;
	//devide the given polynom with the solution we already found.
	root1[1] = first_root;
	poly_after_div[0] = first_poly[0];
	reminder = first_poly[0] * (-root1[1]);
	poly_after_div[1] = first_poly[1] - reminder;
	reminder = poly_after_div[1] * (-root1[1]);
	poly_after_div[2] = first_poly[2] - reminder;
	reminder = poly_after_div[2] * (-root1[1]);
	poly_after_div[3] = first_poly[3] - reminder;
	reminder = poly_after_div[3] * (-root1[1]);
	poly_after_div[4] = first_poly[4] - reminder;
	reminder = poly_after_div[4] * (-root1[1]);

	printf("\nThe new polynom after devide is: ");
	for (i = 0; i < 4; i++)
		printf("%.2lfX^%d+", poly_after_div[i], (4 - i));
	printf("%.2lf\n", poly_after_div[4]);
	
	a1 = poly_after_div[0];
	b1 = poly_after_div[1];
	c1 = poly_after_div[2];
	d1 = poly_after_div[3];
	e1 = poly_after_div[4];

	flag = solve_quartic(a1, b1, c1, d1, e1, &x1, &y1, &x2, &y2, &x3, &y3, &x4, &y4); //this function solve 4th degree polynomial equations.

	//print the solutions of the polynom.
	printf("Equation:%lf x**5 + %lf x**4 + %lf x**3 + %lf x**2 + %lf x +  %lf\n", first_poly[0], first_poly[1], first_poly[2], first_poly[3], first_poly[4], first_poly[5]);
	printf("Solutions are:\n");
	printf("#1) %lf + %lf i\n", first_root, zero);
	printf("#2) %lf + %lf i\n", x1, y1);
	printf("#3) %lf + %lf i\n", x2, y2);
	printf("#4) %lf + %lf i\n", x3, y3);
	printf("#5) %lf + %lf i\n", x4, y4);
	printf("flag = %d\n", flag);

	//check the solutions.
	test_solution(a1, b1, c1, d1, e1, x1, y1, &result1, &result2);
	printf("test for (%lf, %lf) = (%lf, %lf)\n", x1, y1, result1, result2);

	test_solution(a1, b1, c1, d1, e1, x2, y2, &result1, &result2);
	printf("test for (%lf, %lf) = (%lf, %lf)\n", x2, y2, result1, result2);

	test_solution(a1, b1, c1, d1, e1, x3, y3, &result1, &result2);
	printf("test for (%lf, %lf) = (%lf, %lf)\n", x3, y3, result1, result2);

	test_solution(a1, b1, c1, d1, e1, x4, y4, &result1, &result2);
	printf("test for (%lf, %lf) = (%lf, %lf)\n", x4, y4, result1, result2);

	printf("\n\n");

	// This is the second part of the program, we get 6th degree polynom , solve it and calculate the derivative and the minimum points.

	printf("please enter 7 Coefficients,first must be greater then zero:\n ");
	scanf("%lf%lf%lf%lf%lf%lf%lf", &A, &B, &C, &D, &E, &F1, &G);

	while (A <= 0) // if the first coefficient isn't greater then zero.
	{
		printf("Error, the first coefficient must be greater then zero.\n");
		printf("Please enter 7 Coefficients of the polynom (1st coefficient must be greater then zero):\n");
		scanf("%lf%lf%lf%lf%lf%lf%lf", &A, &B, &C, &D, &E, &F1, &G);
	}
	//calculate the derivative of the polynom.
	AX5 = (6 * A);
	BX4 = (5 * B);
	CX3 = (4 * C);
	DX2 = (3 * D);
	EX = (2 * E);
	F = F1;
	//  print the derivative of the polynom.
	printf("the derivative is : ");
	printf("%.3lfX^5+%.3lfX^4+%.3lfX^3+%.3lfX^2+%.3lfX+%.3lf\n\n", AX5, BX4, CX3, DX2, EX, F);

	first_root=BisectionMethod();
	printf("\nBisection returned %.5lf\n", first_root);

	// put the coefficients into another array.
	first_poly[0] = AX5;
	first_poly[1] = BX4;
	first_poly[2] = CX3;
	first_poly[3] = DX2;
	first_poly[4] = EX;
	first_poly[5] = F;

	//devide the given polynom with the solution we already found.
	root1[1] = first_root;
	poly_after_div[0] = first_poly[0];
	reminder = first_poly[0] * (-root1[1]);
	poly_after_div[1] = first_poly[1] - reminder;
	reminder = poly_after_div[1] * (-root1[1]);
	poly_after_div[2] = first_poly[2] - reminder;
	reminder = poly_after_div[2] * (-root1[1]);
	poly_after_div[3] = first_poly[3] - reminder;
	reminder = poly_after_div[3] * (-root1[1]);
	poly_after_div[4] = first_poly[4] - reminder;
	reminder = poly_after_div[4] * (-root1[1]);

	printf("\nThe new polynom after devide is: ");
	for (i = 0; i < 4; i++)
		printf("%.2lfX^%d+", poly_after_div[i], (4 - i));
	printf("%.2lf\n", poly_after_div[4]);
	a1 = poly_after_div[0];
	b1 = poly_after_div[1];
	c1 = poly_after_div[2];
	d1 = poly_after_div[3];
	e1 = poly_after_div[4];

	flag = solve_quartic(a1, b1, c1, d1, e1, &x1, &y1, &x2, &y2, &x3, &y3, &x4, &y4); //this function solve 4th degree polynomial equations.

    //print the solutions of the polynom.
	printf("Equation:%lf x**5 + %lf x**4 + %lf x**3 + %lf x**2 + %lf x +  %lf\n", first_poly[0], first_poly[1], first_poly[2], first_poly[3], first_poly[4], first_poly[5]);
	printf("Solutions are:\n");
	printf("#1) %lf + %lf i\n", first_root, zero);
	printf("#2) %lf + %lf i\n", x1, y1);
	printf("#3) %lf + %lf i\n", x2, y2);
	printf("#4) %lf + %lf i\n", x3, y3);
	printf("#5) %lf + %lf i\n", x4, y4);
	printf("flag = %d\n", flag);
	
	//check the solutions.
	test_solution(a1, b1, c1, d1, e1, x1, y1, &result1, &result2);
	printf("test for (%lf, %lf) = (%lf, %lf)\n", x1, y1, result1, result2);

	test_solution(a1, b1, c1, d1, e1, x2, y2, &result1, &result2);
	printf("test for (%lf, %lf) = (%lf, %lf)\n", x2, y2, result1, result2);

	test_solution(a1, b1, c1, d1, e1, x3, y3, &result1, &result2);
	printf("test for (%lf, %lf) = (%lf, %lf)\n", x3, y3, result1, result2);

	test_solution(a1, b1, c1, d1, e1, x4, y4, &result1, &result2);
	printf("test for (%lf, %lf) = (%lf, %lf)\n", x4, y4, result1, result2);

	// print the minimum value of each solution.
	printf("The minimum points are:\n");
	temp1 = SixthPolynom(first_root);
	printf("p6(%0.4lf) = %.4lf\n",first_root, temp1);
	temp2 = SixthPolynom(x1);
	printf("p6(%0.4lf) = %.4lf\n", x1, temp2);
	temp3 = SixthPolynom(x2);
	printf("p6(%0.4lf) = %.4lf\n", x2, temp3);
	temp4 = SixthPolynom(x3);
	printf("p6(%0.4lf) = %.4lf\n", x3, temp4);
	temp5 = SixthPolynom(x4);
	printf("p6(%0.4lf) = %.4lf\n", x4, temp5);
	
	// calculate the global minimum points and print it.
	min = temp1;
	if (temp2 < min)
		min = temp2;
	if (temp3 <min)
		min = temp3;
	if (temp4 < min)
		min = temp4;
	if (temp2 < min)
		min = temp5;
	printf("Global Minimum value = %0.4lf\n", min);

	system("pause");
	return 0;
}


double SixthPolynom(double x)   //6th degree polynomial function
{
	return A*x*x*x*x*x*x + B*x*x*x*x*x + C*x*x*x*x + D*x*x*x + E*x*x + F*x + G;
}
double FifthPolynom(double x)   //5th degree polynomial function
{
	return AX5 * x*x*x*x*x + BX4 * x*x*x*x + CX3 * x*x*x + DX2 * x*x +EX * x + F;
}

double BisectionMethod() // find the real root with bisection method.
{
	double c1 = -10000, c2 = 10000, solution, d1, d2, d3;

	while (fabsl(c2 - c1) >= EPSILON)
	{
		solution = (c1 + c2) / 2.0;
		d1 = FifthPolynom(c1)*FifthPolynom(solution);
		if (d1 > EPSILON)
			c1 = solution;
		else if (d1 < EPSILON)
			c2 = solution;
	}
	return solution;
}

int solve_quadratic(double a, double b,
	double c,
	double *r1, double *i1,
	double *r2, double *i2)
{
	double delta, delta_root, xv, qv;

	if (a == 0)
		if (b == 0)
			return 0;
		else
		{
			*r1 = -c / b;
			*i1 = 0;
			return 1;
		} /* else */

	xv = -b / (2 * a); /* Take care of square quad */
	qv = a * xv*xv + b * xv + c;
	if (fabs(qv) < EPSILON)
	{
		*r1 = *r2 = xv;
		*i1 = *i2 = 0.0;
		return 2;
	} /* if */


	delta = b * b - 4 * a*c;
	if (delta < 0.0)
	{
		*r1 = *r2 = -b / (2 * a);  /* real part only */
		delta_root = sqrt(-delta);
		*i1 = -delta_root / (2.0*a);
		*i2 = delta_root / (2.0*a);
		return 0;
	} /* if */
	else
	{
		delta_root = sqrt(delta);
		*r1 = (-b - delta_root) / (2.0*a);
		*r2 = (-b + delta_root) / (2.0*a);
		*i1 = 0;
		*i2 = 0;
		return 2;
	} /* else */

} /* solve_quadratic */


void divide_cubic(double a, double *b,
	double *c, double x)
{
	*b = *b + a * x;
	*c = *c + (*b)*x;

} /*  divide_cubic */

void cubic_root(double x, double y, double *r, double *i)
{
	double radius, theta, temp;

	if (x == 0.0 && y == 0.0)
	{
		*r = *i = 0.0;
		return;
	} /* if */


	radius = sqrt(x*x + y * y);
	theta = (asin(fabs(y / radius)));

	temp = exp(log(radius) / 3.0);
	if ((x < 0) && (y >= 0))
		theta = 3 * PI - theta;
	else
		if ((x<0) && (y <= 0))
			theta = 3 * PI + theta;
		else
			if ((x >= 0) && (y<0))
				theta = 4 * PI - theta;

	*r = temp * cos(theta / 3.0);
	*i = temp * sin(theta / 3.0);

} /* cubic_root */

void square_root(double x, double y, double *r, double *i)
{
	double radius, theta, temp;

	if (x == 0.0 && y == 0.0)
	{
		*r = *i = 0.0;
		return;
	} /* if */


	radius = sqrt(x*x + y * y);
	theta = (asin(fabs(y / radius)));

	temp = exp(log(radius) / 2.0);
	if ((x < 0) && (y >= 0))
		theta = 3 * PI - theta;
	else
		if ((x<0) && (y <= 0))
			theta = 3 * PI + theta;
		else
			if ((x >= 0) && (y<0))
				theta = 4 * PI - theta;

	*r = temp * cos(theta / 2.0);
	*i = temp * sin(theta / 2.0);

} /* cubic_root */


int solve_cubic(double a1, double a2,
	double a3, double a4,
	double *x1, double *y1,
	double *x2, double *y2,
	double *x3, double *y3)
{

	double p, q;
	double a, b, c;
	double u1, u2, wr, wi;
	double x, t, ur, ui, vr, vi;
	double b2, c2;
	double i1, i2;
	int n;

	if (a1 == 0.0)
	{
		n = solve_quadratic(a2, a3, a4, x1, &i1, x2, &i2);
		return n;
	} /* if */

	a = a2 / a1;
	b = a3 / a1;
	c = a4 / a1;

	p = b - a * a / 3.0;
	q = c + (2.0*a*a*a - 9.0*a*b) / 27.0;


	if ((p == 0.0) && (q == 0.0))
	{
		*x1 = *x2 = *x3 = -a / 3;
		*y1 = *y2 = *y3 = 0;
		return 3;
	} /* if */

	if (p == 0.0) /* q != 0 */
	{
		cubic_root(-q, 0, &ur, &ui);
		x = ur - (a / 3.0);
		*x1 = x;
		*y1 = 0;

		b2 = a;
		c2 = b;

		divide_cubic(1.0, &b2, &c2, x);
		n = solve_quadratic(1.0, b2, c2, x2, y2, x3, y3);

		return n + 1;

	} /* if */

	solve_quadratic(1.0, -q, -((p*p*p) / 27.0), &u1, &i1, &u2, &i2);

	if (u1 > u2)
	{
		wr = u1;
		wi = i1;
	} /* if */
	else
	{
		wr = u2;
		wi = i2;
	} /* else */

	cubic_root(wr, wi, &ur, &ui);
	cubic_root(-q + wr, wi, &vr, &vi);
	t = vr - ur;
	x = t - (a / 3.0);

	*x1 = x;
	*y1 = 0;

	b2 = a;
	c2 = b;

	divide_cubic(1.0, &b2, &c2, x);
	n = solve_quadratic(1.0, b2, c2, x2, y2, x3, y3);

	return n + 1;

} /* solve_cubic */


void solve_biquadratic(double a, double b, double c,
	double *x1, double *y1,
	double *x2, double *y2,
	double *x3, double *y3,
	double *x4, double *y4)
{

	double q11, q12, q21, q22;
	double sqrt_q11, sqrt_q12, sqrt_q21, sqrt_q22;
	double sqrt_q31, sqrt_q32, sqrt_q41, sqrt_q42;


	solve_quadratic(a, b, c,
		&q11, &q12,
		&q21, &q22);

	square_root(q11, q12, &sqrt_q11, &sqrt_q12);

	sqrt_q21 = -sqrt_q11;
	sqrt_q22 = -sqrt_q12;

	square_root(q21, q22, &sqrt_q31, &sqrt_q32);

	sqrt_q41 = -sqrt_q31;
	sqrt_q42 = -sqrt_q32;

	*x1 = sqrt_q11;
	*y1 = sqrt_q12;

	*x2 = sqrt_q21;
	*y2 = sqrt_q22;


	*x3 = sqrt_q31;
	*y3 = sqrt_q32;


	*x4 = sqrt_q41;
	*y4 = sqrt_q42;


} // solve_biquadratic

int solve_quartic(double a4, double a3,
	double a2, double a1, double a0,
	double *x1, double *y1,
	double *x2, double *y2,
	double *x3, double *y3,
	double *x4, double *y4)
{
	double b, c, d, e;
	double p, q, r;
	double b3, b2, b1, b0;
	double alpha1, alpha2, beta1, beta2, gamma1, gamma2;
	double abs_sqrt_alpha, abs_sqrt_beta;
	double sqrt_alpha1, sqrt_alpha2, sqrt_beta1, sqrt_beta2,
		sqrt_gamma1, sqrt_gamma2;
	double r11, r12, r21, r22, r31, r32, r41, r42;
	double maybe_q;

	int n;

	printf("a4 = %lf, a3 = %lf, a2 = %lf, a1 = %lf, a0 = %lf\n",
		a4, a3, a2, a1, a0);


	b = a3 / a4;
	c = a2 / a4;
	d = a1 / a4;
	e = a0 / a4;

	printf("b = %lf, c = %lf, d = %lf, e = %lf\n", b, c, d, e);


	p = (8.0*c - 3.0*b*b) / 8.0;
	q = (b*b*b - 4.0*b*c + 8.0*d) / 8.0;
	r = (-3.0*b*b*b*b + 256.0*e - 64.0*b*d + 16 * b*b*c) / 256.0;

	printf("p = %lf, q = %lf, r = %lf\n", p, q, r);

	if (fabs(q) < EPSILON)
	{
		solve_biquadratic(1.0, p, r,
			&r11, &r12,
			&r21, &r22,
			&r31, &r32,
			&r41, &r42);

		*x1 = r11 - (b / 4.0);
		*y1 = r12;
		*x2 = r21 - (b / 4.0);
		*y2 = r22;
		*x3 = r31 - (b / 4.0);
		*y3 = r32;
		*x4 = r41 - (b / 4.0);
		*y4 = r42;

		return 1;
	}// if



	b3 = 1.0;
	b2 = 2.0*p;
	b1 = p * p - 4 * r;
	b0 = -q * q;

	printf("b3 = %lf, b2 = %lf, b1 = %lf, b0 = %lf\n", b3, b2, b1, b0);


	n = solve_cubic(b3, b2, b1, b0,
		&alpha1, &alpha2,
		&beta1, &beta2,
		&gamma1, &gamma2);


	printf("alpha1 = %lf, alpha2 = %lf\n", alpha1, alpha2);
	printf("beta1 = %lf, beta2 = %lf\n", beta1, beta2);
	printf("gamma1 = %lf, gamma2 = %lf\n", gamma1, gamma2);

	square_root(alpha1, alpha2, &sqrt_alpha1, &sqrt_alpha2);
	square_root(beta1, beta2, &sqrt_beta1, &sqrt_beta2);
	square_root(gamma1, gamma2, &sqrt_gamma1, &sqrt_gamma2);


	abs_sqrt_alpha = sqrt_alpha1 * sqrt_alpha1 + sqrt_alpha2 * sqrt_alpha2;
	abs_sqrt_beta = sqrt_beta1 * sqrt_beta1 + sqrt_beta2 * sqrt_beta2;


	sqrt_gamma1 = -q * (sqrt_alpha1*sqrt_beta1 -
		sqrt_alpha2 * sqrt_beta2) / (abs_sqrt_alpha*abs_sqrt_beta);
	sqrt_gamma2 = q * (sqrt_alpha2*sqrt_beta1 +
		sqrt_alpha1 * sqrt_beta2) / (abs_sqrt_alpha*abs_sqrt_beta);


	printf("sqrt_alpha1 = %lf, sqrt_alpha2 = %lf\n", sqrt_alpha1, sqrt_alpha2);
	printf("sqrt_beta1 = %lf, sqrt_beta2 = %lf\n", sqrt_beta1, sqrt_beta2);
	printf("sqrt_gamma1 = %lf, sqrt_gamma2 = %lf\n", sqrt_gamma1, sqrt_gamma2);

	r11 = (sqrt_alpha1 + sqrt_beta1 + sqrt_gamma1) / 2.0;
	r12 = (sqrt_alpha2 + sqrt_beta2 + sqrt_gamma2) / 2.0;
	r21 = (sqrt_alpha1 - sqrt_beta1 - sqrt_gamma1) / 2.0;
	r22 = (sqrt_alpha2 - sqrt_beta2 - sqrt_gamma2) / 2.0;
	r31 = (-sqrt_alpha1 + sqrt_beta1 - sqrt_gamma1) / 2.0;
	r32 = (-sqrt_alpha2 + sqrt_beta2 - sqrt_gamma2) / 2.0;
	r41 = (-sqrt_alpha1 - sqrt_beta1 + sqrt_gamma1) / 2.0;
	r42 = (-sqrt_alpha2 - sqrt_beta2 + sqrt_gamma2) / 2.0;

	printf("r11 = %lf, r12 = %lf, r21 = %lf,r22 = %lf\n",
		r11, r12, r21, r22);
	printf("r31 = %lf, r32 = %lf, r41 = %lf,r42 = %lf\n",
		r31, r32, r41, r42);


	*x1 = r11 - (b / 4.0);
	*y1 = r12;
	*x2 = r21 - (b / 4.0);
	*y2 = r22;
	*x3 = r31 - (b / 4.0);
	*y3 = r32;
	*x4 = r41 - (b / 4.0);
	*y4 = r42;

	printf("b/4.0 = %lf\n", b / 4.0);

	return 2;
} // solve_quartic

void complex_mult(double a, double b,
	double c, double d, double *x, double *y)
{
	*x = a * c - b * d;
	*y = a * d + c * b;
} // complex_mult

void complex_add(double a, double b,
	double c, double d, double *x, double *y)
{
	*x = a + b;
	*y = c + d;
} // complex_mult

void test_solution(double a, double b,
	double c, double d, double e, double x, double y,
	double *result1, double *result2)
{
	double x1, y1, temp11, temp12, temp21, temp22;
	int i;

	temp11 = x;
	temp12 = y;
	temp21 = e;
	temp22 = 0.0;

	temp21 += d * temp11;
	temp22 += d * temp12;;

	complex_mult(temp11, temp12, x, y, &temp11, &temp12);
	temp21 += c * temp11;
	temp22 += c * temp12;;

	complex_mult(temp11, temp12, x, y, &temp11, &temp12);
	temp21 += b * temp11;
	temp22 += b * temp12;;

	complex_mult(temp11, temp12, x, y, &temp11, &temp12);
	temp21 += a * temp11;
	temp22 += a * temp12;;

	*result1 = temp21;
	*result2 = temp22;

} // test_solution 