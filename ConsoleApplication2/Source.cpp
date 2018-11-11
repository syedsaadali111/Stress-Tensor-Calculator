#include <stdio.h>
#include <math.h>
#include <stdlib.h>


typedef struct { double r, a, b, c; } tEquation;

void calculateRoots(double a, double b, double c, double d, double *x1, double *x2, double *x3)
{
	double e, f, g, h, i, j, k, l, m, n, p, r, s, t, u;
	int w;
	e = 2.7182818284590;
	f = ((3.0 * c / a) - (b*b / (a*a))) / 3.0;   // ** bracketed (a*a)!
	g = (((2.0 * b*b*b) / (a*a*a)) - (9.0 * b * c / (a*a)) + (27.0 * d / a)) / 27.0;   // ** brackets!
	h = (g*g / 4.0) + (f*f*f / 27.0);
	i = sqrtf(((g*g / 4) - h));
	j = expf(log10(i) / log10(e) / 3);
	k = acos((-1)*(g / (2 * i)));
	l = j*(-1);
	m = cos(k / 3);
	n = sqrtf(3)*sin(k / 3.0);
	p = (b / (3.0 * a))*(-1);
	r = (-1)*(g / 2) + sqrtf(h);
	s = expf(log10(r) / log10(e) / 3);
	t = (-1)*(g / 2) - sqrt(h);
	u = expf(log10(t) / log10(e) / 3);
	if (h>0) w = 1;
	if (h <= 0) w = 3;
	if ((f == 0) && (g == 0) && (h == 0)) w = 2;
	switch (w) {
	case 1:
		*x1 = (s + u) - (b / 3.0 * a);
		*x2 = (-1)*(s + u) / 2 - (b / 3 * a);
		*x3 = (s - u)*sqrt(3) / 2;
		break;
	case 2:
		*x1 = exp(log10(d / a) / log10(e) / 3)*(-1);
		break;
	case 3:
		*x1 = 2 * j*cos(k / 3) - (b / (3 * a));
		*x2 = l*(m + n) + p;
		*x3 = l*(m - n) + p;
		break;
	}
}

void simultaneousEqSol(tEquation equ1[], double *a, double *b, double *c)
{
	tEquation equ2[2], equ3[1];

	equ2[0].r = equ1[0].r * equ1[1].c - equ1[1].r * equ1[0].c;
	equ2[0].a = equ1[0].a * equ1[1].c - equ1[1].a * equ1[0].c;
	equ2[0].b = equ1[0].b * equ1[1].c - equ1[1].b * equ1[0].c;
	equ2[0].c = 0;

	equ2[1].r = equ1[1].r * equ1[2].c - equ1[2].r * equ1[1].c;
	equ2[1].a = equ1[1].a * equ1[2].c - equ1[2].a * equ1[1].c;
	equ2[1].b = equ1[1].b * equ1[2].c - equ1[2].b * equ1[1].c;
	equ2[1].c = 0;

	equ3[0].r = equ2[0].r * equ2[1].b - equ2[1].r * equ2[0].b;
	equ3[0].a = equ2[0].a * equ2[1].b - equ2[1].a * equ2[0].b;
	equ3[0].b = 0;
	equ3[0].c = 0;

	*a = equ3[0].r / equ3[0].a;

	*b = (equ2[0].r - equ2[0].a * *a) / equ2[0].b;

	*c = (equ1[0].r - equ1[0].a * *a - equ1[0].b * *b) / equ1[0].c;

}

void l_m_nCalculator(double x, double y, double z, double *l, double *m, double *n)
{

	double vectorMagnitude = sqrtf(powf(x, 2) + powf(y, 2) + powf(z, 2));
	*l = x / vectorMagnitude * 1.0;
	*m = y / vectorMagnitude * 1.0;
	*n = z / vectorMagnitude * 1.0;
}

void swap(double *a, double *b)
{
	double temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

void sortX(double *x1, double *x2, double *x3)
{
	if (*x1 < *x3)
		swap(x1, x3);
	if (*x1 < *x2)
		swap(x1, x2);
	if (*x2 < *x3)
		swap(x2, x3);
}

int main()
{
	double x1, x2, x3, a, b, c, d,
		Sx,
		Sy,
		Sz,
		Sxy,
		Sxz,
		Syz,
		i1, i2, i3;
	double l, m, n;
	double ax, ay, az;
	
	printf("Stress Tensor:\n\n");
	printf("\nEnter normal x:");
	scanf("%lf", &Sx);
	printf("\nEnter normal y:");
	scanf("%lf", &Sy);
	printf("\nEnter normal z:");
	scanf("%lf", &Sz);
	printf("\nEnter shear xy:");
	scanf("%lf", &Sxy);
	printf("\nEnter shear xz:");
	scanf("%lf", &Sxz);
	printf("\nEnter shear yz:");
	scanf("%lf", &Syz);

	i1 = Sx + Sy + Sz;
	i2 = (Sx * Sy + Sy * Sz + Sz * Sx) - (powf(Sxy, 2) + powf(Syz, 2) + powf(Sxz, 2));
	i3 = (Sx * Sy * Sz) + 2 * Sxy*Syz*Sxz - (Sx * powf(Syz, 2) + Sy*powf(Sxz, 2) + Sz * powf(Sxy, 2));

	//printf("\n\ni1: %f, i2: %f, i3: %f", i1, i2, i3);

	calculateRoots(1, -i1, i2, -i3, &x1, &x2, &x3);

	sortX(&x1, &x2, &x3);
	

	/*x1 = ((int)(x1 * 10)) / 10.0;
	x2 = ((int)(x2 * 10)) / 10.0;
	x3 = ((int)(x3 * 10)) / 10.0;*/

	printf("\n\nx1: %f\nx2: %f\nx3: %f\n", x1, x2, x3);

	// **********FOR X1********** 

	/* tEquation equ1[] = {
		{ 0,  (Sx - x1), Sxy, Sxz },
		{ 0, Sxy, (Sy - x1), Syz },
		{ 1, 0, 0, 1},
	}; */

	//simultaneousEqSol(equ1, &ax, &ay, &az);

	printf("\nX1:\n");

	ax = ((Sy - x1) * (Sz - x1)) - (Syz * Syz); 
	ay = -(Sxy * (Sz - x1) - (Sxz * Syz));
	az = (Sxy * Syz) - (Sxz * (Sy - x1));

	//printf("Ax: %f\nAy: %f\nAz: %f", ax,ay,az);

	l_m_nCalculator(ax, ay, az, &l, &m, &n);

	printf("\nL: %f\nM: %f\nN: %f\n", l, m, n);

	//**************FOR X2**************

	printf("\nX2:\n");

	ax = -(Sxy * (Sz - x2) - (Sxz * Syz));
	ay = ((Sx - x2) * (Sz - x2)) - (Sxz *Sxz);
	az = -((Sx - x2) * Syz - (Sxz * Sxy));

	//printf("\nAx: %f\nAy: %f\nAz: %f", ax, ay, az);

	l_m_nCalculator(ax, ay, az, &l, &m, &n);

	printf("\nL: %f\nM: %f\nN: %f\n", l, m, n);

	//**************FOR X3*******************
	
	printf("\nX3:\n");

	ax = Sxy * Syz - (Sy - x3) * Sxz;
	ay = -((Sx - x3) * Syz - (Sxy * Sxz));
	az = (Sx - x3) * (Sy - x3) - (Sxy * Sxy);

	//printf("\nAx: %f\nAy: %f\nAz: %f", ax, ay, az);

	l_m_nCalculator(ax, ay, az, &l, &m, &n);

	printf("\nL: %f\nM: %f\nN: %f\n", l, m, n);

	system("PAUSE");
	return 0;


}
