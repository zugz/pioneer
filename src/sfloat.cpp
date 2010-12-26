#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include "sfloat.h"

sfloat sfloat::Log(sfloat a)
{
	const sfloat ln2(2977044471U,-32,false);
	if (a <= sfloat(0,0,false)) {
		return sfloat(0,0,false,TYPE_NAN);
	}
	// do powers of two
	sfloat r(0,0,false);
	while (a > sfloat(1,0,false)) {
		r += ln2; // log(2)
		a *= sfloat(1,-1,false);
	}
	while (a < sfloat(1,-1,false)) {
		r -= ln2; // log(0.5)
		a *= sfloat(2,0,false);
	}
	// do remaining fraction
	sfloat x = sfloat(1,0,false) - a;
	sfloat x0 = x;
	r -= x;
	for (int i=2; i<28; i++) {
		x = x*x0;
		r -= x * inverses[i];
	}
	return r;
}

sfloat sfloat::Exp(sfloat val)
{
	sfloat fraction = val;
	sfloat rhigh = sfloat(1,0,false);
	// do e^(+ve power of two) bits
	if (fraction > sfloat(1,0,false)) {
		sfloat intpart = fraction.floor();
		if (intpart.m >= 1024) {
			// would overflow 16-bit exponent
			return sfloat(0,0,false,TYPE_INF);
		}
		fraction -= intpart;
		for(int i=9; i>=0; i--) {
			if (intpart.m & (1<<i)) {
				rhigh *= exp_pow2[i];
			}
		}
	}
	// do e^(+ve power of two) bits
	if (fraction < sfloat(1,0,true)) {
		sfloat intpart = fraction.floor();
		fraction -= intpart;
		for(int i=9; i>=0; i--) {
			if (intpart.m & (1<<i)) {
				rhigh *= exp_negpow2[i];
			}
		}
	}
	sfloat r(1,0,false);
	sfloat x = fraction;
	r += fraction;
	for (int i=2; i<12; i++) {
		x = x * fraction;
		r += x * inverseFactorials[i];
	}
	return r * rhigh;
}

const sfloat sfloat::oneThird(2863311530U,-33,false);
const sfloat sfloat::inverseFactorials[32] = {
	sfloat(2147483648U,-31,false),
	sfloat(2147483648U,-31,false),
	sfloat(2147483648U,-32,false),
	sfloat(2863311530U,-34,false),
	sfloat(2863311530U,-36,false),
	sfloat(2290649224U,-38,false),
	sfloat(3054198966U,-41,false),
	sfloat(3490513104U,-44,false),
	sfloat(3490513104U,-47,false),
	sfloat(3102678314U,-50,false),
	sfloat(2482142651U,-53,false),
	sfloat(3610389311U,-57,false),
	sfloat(2406926207U,-60,false),
	sfloat(2962370717U,-64,false),
	sfloat(3385566534U,-68,false),
	sfloat(3611270969U,-72,false),
	sfloat(3611270969U,-76,false),
	sfloat(3398843266U,-80,false),
	sfloat(3021194015U,-84,false),
	sfloat(2544163381U,-88,false),
	sfloat(4070661411U,-93,false),
	sfloat(3101456313U,-97,false),
	sfloat(2255604591U,-101,false),
	sfloat(3138232475U,-106,false),
	sfloat(4184309968U,-111,false),
	sfloat(2677958380U,-115,false),
	sfloat(3295948776U,-120,false),
	sfloat(3906309661U,-125,false),
	sfloat(2232176949U,-129,false),
	sfloat(2463091806U,-134,false),
	sfloat(2627297927U,-139,false)
};
const sfloat sfloat::inverses[32] = {
	sfloat(2147483648U,-31,false),
	sfloat(2147483648U,-31,false),
	sfloat(2147483648U,-32,false),
	sfloat(2863311530U,-33,false),
	sfloat(2147483648U,-33,false),
	sfloat(3435973836U,-34,false),
	sfloat(2863311530U,-34,false),
	sfloat(2454267026U,-34,false),
	sfloat(2147483648U,-34,false),
	sfloat(3817748707U,-35,false),
	sfloat(3435973836U,-35,false),
	sfloat(3123612578U,-35,false),
	sfloat(2863311530U,-35,false),
	sfloat(2643056797U,-35,false),
	sfloat(2454267026U,-35,false),
	sfloat(2290649224U,-35,false),
	sfloat(2147483648U,-35,false),
	sfloat(4042322160U,-36,false),
	sfloat(3817748707U,-36,false),
	sfloat(3616814565U,-36,false),
	sfloat(3435973836U,-36,false),
	sfloat(3272356035U,-36,false),
	sfloat(3123612578U,-36,false),
	sfloat(2987803336U,-36,false),
	sfloat(2863311530U,-36,false),
	sfloat(2748779069U,-36,false),
	sfloat(2643056797U,-36,false),
	sfloat(2545165805U,-36,false),
	sfloat(2454267026U,-36,false),
	sfloat(2369637128U,-36,false),
	sfloat(2290649224U,-36,false)
};
const sfloat sfloat::exp_pow2[11] = {
	sfloat(2918732888U,-30,false),
	sfloat(3966969286U,-29,false),
	sfloat(3664019825U,-26,false),
	sfloat(3125761002U,-20,false),
	sfloat(2274844293U,-8,false),
	sfloat(2409758306U,15,false),
	sfloat(2704064871U,61,false),
	sfloat(3404899886U,153,false),
	sfloat(2699285568U,338,false),
	sfloat(3392874533U,707,false),
	sfloat(2147483648U,993,false)
};
const sfloat sfloat::exp_negpow2[10] = {
	sfloat(3160060337U,-33,false), // exp(0)
	sfloat(2325042461U,-34,false), // exp(-1)
	sfloat(2517282241U,-37,false),
	sfloat(2950760480U,-43,false),
	sfloat(4054506967U,-55,false),
	sfloat(3827509179U,-78,false),
	sfloat(3410928537U,-124,false),
	sfloat(2708852636U,-216,false),
	sfloat(3416967861U,-401,false),
	sfloat(2718453613U,-770,false),
};

#ifdef TEST
static void testLog()
{
	printf("Testing log(): ");
	double min = FLT_MAX;
	double max = -FLT_MAX;
	for (int i=-16; i<10000; i++) {
		sfloat v(rand()&0x7fffffff, (rand()%500) - 272, false);
		if (v.ToDouble() < min) min = v.ToDouble();
		if (v.ToDouble() > max) max = v.ToDouble();
		double model = log(v.ToDouble());
		double test = sfloat::Log(v).ToDouble();
	//	printf("Log %e = %e (%e)\n", v.ToDouble(), test, model);
		assert ((model == test) || (fabs(test/model - 1.0) < 0.0000001));
	}
	printf("ok (range %e - %e)\n", min, max);
}
static void testExp()
{
	printf("Testing exp(): ");
	double min = FLT_MAX;
	double max = -FLT_MAX;
	for (int i=1; i<10000; i++) {
		sfloat v;
		if (rand()&1) {
			v = sfloat(rand()%(1<<19), -10, true);
		} else {
			v = sfloat(rand()%(1024*1024), -10, false);
		}
		double model = exp(v.ToDouble());
		double test = sfloat::Exp(v).ToDouble();
		if (v.ToDouble() < min) min = v.ToDouble();
		if (v.ToDouble() > max) max = v.ToDouble();
	//	printf("exp %e = %e (%e)\n", v.ToDouble(), test, model);
		assert ((model == test) || (fabs(test/model - 1.0) < 0.00000001));
	}
	printf("ok (range %e - %e)\n", min, max);
}
static void testMul()
{
	printf("Testing mul(): ");
	for (int i=1; i<2000; i++) {
		sfloat a(rand(),(rand()&0x1f)-0xf,rand()&1 ? true : false);
		sfloat b(rand(),(rand()&0x1f)-0xf,rand()&1 ? true : false);
		if ((rand()&0x1f) == 0) a = sfloat(0,0,false);
		if ((rand()&0x1f) == 0) b = sfloat(0,0,false);
		double model = a.ToDouble() * b.ToDouble();
		double test = (a*b).ToDouble();
		//printf("mul %e*%e = %e (%e)\n", a.ToDouble(), b.ToDouble(), test, model);
		assert((a*b)==(b*a));
		assert ((model == test) || fabs(test/model - 1.0) < 0.000000001);
	}
	printf("ok\n");
}
static void testDiv()
{
	printf("Testing div(): ");
	for (int i=1; i<2000; i++) {
		sfloat a(rand(),(rand()&0x1f)-48,rand()&1 ? true : false);
		sfloat b(rand(),(rand()&0x1f)-48,rand()&1 ? true : false);
		if ((rand()&0x1f) == 0) a = sfloat(0,0,false);
		double model = a.ToDouble() / b.ToDouble();
		double test = (a/b).ToDouble();
		//printf("div %e*%e = %e (%e)\n", a.ToDouble(), b.ToDouble(), test, model);
		assert ((model == test) || fabs(test/model - 1.0) < 0.000000001);
	}
	printf("ok\n");
}
static void testAdd()
{
	printf("Testing add(): ");
	for (int i=1; i<2000; i++) {
		sfloat a(rand()<<1,rand()&0xf,rand()&1 ? true : false);
		sfloat b(rand()<<1,rand()&0xf,rand()&1 ? true : false);
		if ((rand()&0x1f) == 0) a = sfloat(0,0,false);
		if ((rand()&0x1f) == 0) b = sfloat(0,0,false);
		double model = a.ToDouble() + b.ToDouble();
		double test = (a+b).ToDouble();
		//printf("add %e+%e = %e (%e)\n", a.ToDouble(), b.ToDouble(), test, model);
		assert((a+b)==(b+a));
		assert ((model == test) || fabs(test/model - 1.0) < 0.000000001);
	}
	printf("ok\n");
}
static void testSub()
{
	printf("Testing sub(): ");
	for (int i=1; i<2000; i++) {
		sfloat a(rand()<<1,rand()&0xf,rand()&1 ? true : false);
		sfloat b(rand()<<1,rand()&0xf,rand()&1 ? true : false);
		if ((rand()&0x1f) == 0) a = sfloat(0,0,false);
		if ((rand()&0x1f) == 0) b = sfloat(0,0,false);
		double model = a.ToDouble() - b.ToDouble();
		double test = (a-b).ToDouble();
		//printf("sub %e-%e = %e (%e)\n", a.ToDouble(), b.ToDouble(), test, model);
		assert ((model == test) || fabs(test/model - 1.0) < 0.000000001);
	}
	printf("ok\n");
}
static void testCmp()
{
	printf("Testing cmp(): ");
	for (int i=1; i<20000; i++) {
		sfloat a(rand()<<1,rand()&0xf,rand()&1 ? true : false);
		sfloat b(rand()<<1,rand()&0xf,rand()&1 ? true : false);
		if ((rand()&0x1f) == 0) a = sfloat(0,0,false);
		if ((rand()&0x1f) == 0) b = sfloat(0,0,false);
	//	printf("%f > %f\n", a.ToDouble(), b.ToDouble());
		assert ((a.ToDouble() > b.ToDouble()) == (a > b));
		assert ((a.ToDouble() >= b.ToDouble()) == (a >= b));
		assert ((a.ToDouble() <= b.ToDouble()) == (a <= b));
		assert ((a.ToDouble() < b.ToDouble()) == (a < b));
		assert ((a.ToDouble() == b.ToDouble()) == (a == b));
	}
	printf("ok\n");
}

/* Don't use this shit for anything since it won't
 * give identical results on all platforms (probably) */
static sfloat ieee_double_to_sfloat(double f)
{
	Uint64 data = *((Uint64*)&f);
	Uint64 m = (data & 0xfffffffffffffULL) | (1ULL<<52);
	m >>= 21;
	short e = ((data>>52) & 0x7ff) - 1023 - 31;
	char sign = (data>>63) & 1;
	return sfloat(m, e, sign);
}

int main()
{
	srand(1234);
	testAdd();
	testSub();
	testMul();
	testDiv();
	testCmp();
	testLog();
	testExp();

//	ieee_double_to_sfloat(0.33333333333333333333333333333333333333333333333).Print();
//	ieee_double_to_sfloat(log(2.0)).Print();
//	printf("%f\n", sfloat::Cuberoot(sfloat(36,0,true)).ToDouble());

	return 0;
}
#endif /* TEST */
