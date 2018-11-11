#include <stdio.h>

typedef struct { double r, a, b, c; } tEquation;
tEquation equ1[] = {
	{ 0,  40, 50, 60 },      // -44.3940 = 50a + 37b + c (A)
	{ 0, 60, 40, 50 },      // -45.3049 = 43a + 39b + c (B)
	{ 0,  50, 40, 60 },      // -44.9594 = 52a + 41b + c (C)
};
tEquation equ2[2], equ3[1];

static void dumpEqu(char *desc, tEquation *e, char *post) {
	printf("%10s: %12.8lf = %12.8lfa + %12.8lfb + %12.8lfc (%s)\n",
		desc, e->r, e->a, e->b, e->c, post);
}

int main(void) {
	double a, b, c;
	dumpEqu(">", &(equ1[0]), "A");
	dumpEqu(">", &(equ1[1]), "B");
	dumpEqu(">", &(equ1[2]), "C");
	puts("");

	// A - B
	equ2[0].r = equ1[0].r * equ1[1].c - equ1[1].r * equ1[0].c;
	equ2[0].a = equ1[0].a * equ1[1].c - equ1[1].a * equ1[0].c;
	equ2[0].b = equ1[0].b * equ1[1].c - equ1[1].b * equ1[0].c;
	equ2[0].c = 0;

	// B - C
	equ2[1].r = equ1[1].r * equ1[2].c - equ1[2].r * equ1[1].c;
	equ2[1].a = equ1[1].a * equ1[2].c - equ1[2].a * equ1[1].c;
	equ2[1].b = equ1[1].b * equ1[2].c - equ1[2].b * equ1[1].c;
	equ2[1].c = 0;

	dumpEqu("A-B", &(equ2[0]), "D");
	dumpEqu("B-C", &(equ2[1]), "E");
	puts("");

	equ3[0].r = equ2[0].r * equ2[1].b - equ2[1].r * equ2[0].b;
	equ3[0].a = equ2[0].a * equ2[1].b - equ2[1].a * equ2[0].b;
	equ3[0].b = 0;
	equ3[0].c = 0;

	dumpEqu("D-E", &(equ3[0]), "F");
	puts("");

	a = equ3[0].r / equ3[0].a;
	printf("From (F    ), a = %12.8lf (G)\n", a);

	b = (equ2[0].r - equ2[0].a * a) / equ2[0].b;
	printf("From (D,G  ), b = %12.8lf (H)\n", b);

	c = (equ1[0].r - equ1[0].a * a - equ1[0].b * b) / equ1[0].c;
	printf("From (A,G,H), c = %12.8lf (I)\n", c);

	return 0;
}