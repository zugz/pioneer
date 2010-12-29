#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include "sfloat.h"

#define SFLOAT_ONE  sfloat(0x3f800000)
#define SFLOAT_HALF sfloat(0x3f000000)
#define SFLOAT_NEG_ONE  sfloat(0xbf800000)

sfloat sfloat::Log(sfloat a)
{
	const sfloat ln2(0x3f317218); // log(2)
	if (a <= sfloat(0)) {
		sfloat r;
		r.val = SFLOAT_NVE_INF;
		return r;
	}
	// do powers of two
	sfloat r(0);
	while (a > SFLOAT_ONE) {
		r += ln2; // log(2)
		a.val -= (1<<23); // divide by two
	}
	while (a < SFLOAT_HALF) {
		r -= ln2; // log(0.5)
		a.val += (1<<23); // multiply by two
	}
	// do remaining fraction
	sfloat x = SFLOAT_ONE - a;
	sfloat x0 = x;
	r -= x;
	for (int i=2; i<28; i++) {
		x = x*x0;
		r -= x * inverses[i];
	}
	return r;
}

sfloat sfloat::Exp(sfloat p)
{
	// p > 256
	if (p > sfloat(0x43800000)) {
		sfloat r;
		r.val = SFLOAT_PVE_INF;
		return r;
	}
	if (p < sfloat(0x42b171ab)) return sfloat(0);
	sfloat fraction = p;
	sfloat rhigh = SFLOAT_ONE;
	// do e^(+ve power of two) bits
	if (fraction > SFLOAT_ONE) {
		int intpart = fraction.ToInt32();
		fraction -= sfloat(intpart,1);
		for(int i=8; i>=0; i--) {
			if (intpart & (1<<i)) {
				rhigh *= exp_pow2[i];
			}
		}
	}
	// do -e^(+ve power of two) bits
	if (fraction < SFLOAT_NEG_ONE) {
		int intpart = -fraction.ToInt32();
		fraction += sfloat(intpart,1);
		for(int i=7; i>=0; i--) {
			if (intpart & (1<<i)) {
				rhigh *= exp_negpow2[i];
			}
		}
	}
	sfloat r = SFLOAT_ONE;
	sfloat x = fraction;
	r += fraction;
	for (int i=2; i<11; i++) {
		x = x * fraction;
		r += x * inverseFactorials[i];
	}
	return r * rhigh;
}

const sfloat sfloat::oneThird(2863311530U,-33,false);
const sfloat sfloat::inverseFactorials[32] = {
	sfloat(0x3f800000),
	sfloat(0x3f800000),
	sfloat(0x3f000000),
	sfloat(0x3e2aaaab),
	sfloat(0x3d2aaaab),
	sfloat(0x3c088889),
	sfloat(0x3ab60b61),
	sfloat(0x39500d01),
	sfloat(0x37d00d01),
	sfloat(0x3638ef1d),
	sfloat(0x3493f27e),
	sfloat(0x32d7322b),
	sfloat(0x310f76c7),
	sfloat(0x2f309231),
	sfloat(0x2d49cba6),
	sfloat(0x2b573fa0),
	sfloat(0x29573fa0),
	sfloat(0x274a963c),
	sfloat(0x253413c3),
	sfloat(0x2317a4db),
	sfloat(0x20f2a15d),
	sfloat(0x1eb8dc78),
	sfloat(0x1c8671cb),
	sfloat(0x1a3b0da0),
	sfloat(0x17f9677f),
	sfloat(0x159f9e66),
	sfloat(0x1344742f),
	sfloat(0x10e8d58d),
	sfloat(0x0e850c50),
	sfloat(0x0c12cfcb),
	sfloat(0x099c9961),
	sfloat(0x0721a696),
};
const sfloat sfloat::inverses[32] = {
	sfloat(0x3f800000),
	sfloat(0x3f800000),
	sfloat(0x3f000000),
	sfloat(0x3eaaaaab),
	sfloat(0x3e800000),
	sfloat(0x3e4ccccd),
	sfloat(0x3e2aaaab),
	sfloat(0x3e124925),
	sfloat(0x3e000000),
	sfloat(0x3de38e39),
	sfloat(0x3dcccccd),
	sfloat(0x3dba2e8c),
	sfloat(0x3daaaaab),
	sfloat(0x3d9d89d9),
	sfloat(0x3d924925),
	sfloat(0x3d888889),
	sfloat(0x3d800000),
	sfloat(0x3d70f0f1),
	sfloat(0x3d638e39),
	sfloat(0x3d579436),
	sfloat(0x3d4ccccd),
	sfloat(0x3d430c31),
	sfloat(0x3d3a2e8c),
	sfloat(0x3d321643),
	sfloat(0x3d2aaaab),
	sfloat(0x3d23d70a),
	sfloat(0x3d1d89d9),
	sfloat(0x3d17b426),
	sfloat(0x3d124925),
	sfloat(0x3d0d3dcb),
	sfloat(0x3d088889),
	sfloat(0x3d042108),
};
const sfloat sfloat::exp_pow2[8] = {
	sfloat(0x402df854),
	sfloat(0x40ec7326),
	sfloat(0x425a6481),
	sfloat(0x453a4f54),
	sfloat(0x4b07975f),
	sfloat(0x568fa1fe),
	sfloat(0x6da12cc1),
	sfloat(0x7f800000),
};
const sfloat sfloat::exp_negpow2[7] = {
	sfloat(0x3ebc5ab2), // exp(-1)
	sfloat(0x3e0a9555), // exp(-2)
	sfloat(0x3c960aae), // exp(-4)
	sfloat(0x39afe108),
	sfloat(0x33f1aade),
	sfloat(0x28642328),
	sfloat(0x114b4ea4),
};

#ifdef TEST
static sfloat RandSfloat() {
	Uint32 r;
	do {
		r = (rand()<<16) | (rand()&0xffff);
	} while (((r>>23)&0xff)==0xff); // ignore INF and NaN
	return sfloat(r);
}

static void testLog()
{
	printf("Testing log(): ");
	float min = FLT_MAX;
	float max = -FLT_MAX;
	for (int i=-16; i<1000; i++) {
		sfloat v(rand()&0x7fffffff, (rand()%400) - 222, false);
		if (v.ToFloat() < min) min = v.ToFloat();
		if (v.ToFloat() > max) max = v.ToFloat();
		float model = log(v.ToFloat());
		float test = sfloat::Log(v).ToFloat();
		printf("Log %e = %e (%e)\n", v.ToFloat(), test, model);
		assert ((model == test) || (fabs(test/model - 1.0) < 0.00001));
	}
	printf("ok (range %e - %e)\n", min, max);
}

static void testExp()
{
	printf("Testing exp(): ");
	float min = FLT_MAX;
	float max = -FLT_MAX;
	for (int i=1; i<100000; i++) {
		sfloat v = RandSfloat();/*
		if (rand()&1) {
			v = sfloat(rand()%(1<<19), -15, true);
		} else {
			v = sfloat(rand()%(1024*1024), -10, false);
		}*/
		float model = exp(v.ToFloat());
		float test = sfloat::Exp(v).ToFloat();
		if (v.ToFloat() < min) min = v.ToFloat();
		if (v.ToFloat() > max) max = v.ToFloat();
		printf("exp %e = %e (%e)\n", v.ToFloat(), test, model);
		assert ((model == test) || (fabs(test/model - 1.0) < 0.00001));
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
		float model = a.ToFloat() * b.ToFloat();
		float test = (a*b).ToFloat();
		//printf("mul %e*%e = %e (%e)\n", a.ToFloat(), b.ToFloat(), test, model);
		assert((a*b)==(b*a));
		assert ((model == test) || fabs(test/model - 1.0) < 0.000001);
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
		float model = a.ToFloat() / b.ToFloat();
		float test = (a/b).ToFloat();
		//printf("div %e*%e = %e (%e)\n", a.ToFloat(), b.ToFloat(), test, model);
		assert ((model == test) || fabs(test/model - 1.0) < 0.000001);
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
		float model = a.ToFloat() + b.ToFloat();
		float test = (a+b).ToFloat();
		//printf("add %e+%e = %e (%e)\n", a.ToFloat(), b.ToFloat(), test, model);
		assert((a+b)==(b+a));
		assert ((model == test) || fabs(test/model - 1.0) < 0.000001);
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
		float model = a.ToFloat() - b.ToFloat();
		float test = (a-b).ToFloat();
		//printf("sub %e-%e = %e (%e)\n", a.ToFloat(), b.ToFloat(), test, model);
		assert ((model == test) || fabs(test/model - 1.0) < 0.000001);
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
	//	printf("%f > %f\n", a.ToFloat(), b.ToFloat());
		assert ((a.ToFloat() > b.ToFloat()) == (a > b));
		assert ((a.ToFloat() >= b.ToFloat()) == (a >= b));
		assert ((a.ToFloat() <= b.ToFloat()) == (a <= b));
		assert ((a.ToFloat() < b.ToFloat()) == (a < b));
		assert ((a.ToFloat() == b.ToFloat()) == (a == b));
	}
	printf("ok\n");
}

int main()
{
	srand(1234);
	testAdd();
	testSub();
	testMul();
	testDiv();
	testCmp();
	testExp();
	testLog();

	return 0;
}
#endif /* TEST */
