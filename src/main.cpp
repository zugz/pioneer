#include "libs.h"
#include "Pi.h"
#include <signal.h>

void sigsegv_handler(int signum)
{
	if (signum == SIGSEGV) {
		printf("Segfault! All is lost! Abandon ship!\n");
		SDL_Quit();
		abort();
	}
}
#include "fixed.h"

int main(int argc, char**)
{
	/*
	printf("%.9f, %.9f\n", fixed::Exp(fixed(31,1)).ToDouble(), exp(31.0));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(214,10)).ToDouble(), exp(21.4));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(205,10)).ToDouble(), exp(20.5));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(20,1)).ToDouble(), exp(20.0));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(189,10)).ToDouble(), exp(18.9));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(183,10)).ToDouble(), exp(18.3));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(165,10)).ToDouble(), exp(16.5));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(84,10)).ToDouble(), exp(8.4));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(42,10)).ToDouble(), exp(4.2));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(21,10)).ToDouble(), exp(2.1));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(105,100)).ToDouble(), exp(1.05));
	printf("%.9f, %.9f\n", fixed::Exp(fixed(1,2)).ToDouble(), exp(0.5));
	exit(0);
	*/
//	signal(SIGSEGV, sigsegv_handler);
	Pi::Init();
	for (;;) Pi::Start();
	return 0;
}
