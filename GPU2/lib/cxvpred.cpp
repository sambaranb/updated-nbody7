#include <cstdio>
#include <cassert>

#define CHKALIGN(ptr, num) assert((unsigned long)(ptr) % num == 0);
#define _out_ 
#define inout 
#define RCPTR(x) (* __restrict const x)

void cxvpred(
		const int js,
		const int je,
		const double time,
		const double RCPTR(t0   ),
		const double RCPTR(x0   )[3],
		const double RCPTR(x0dot)[3],
		const double RCPTR(f    )[3],
		const double RCPTR(fdot )[3],
		_out_ double RCPTR(x    )[3],
		_out_ double RCPTR(xdot )[3],
		inout double RCPTR(tpred))
{
	// fprintf(stderr, "CXVPRED : js = %d, je = %d\n", js, je);
	CHKALIGN(js, 2);
	CHKALIGN(je, 2);
#pragma omp parallel for
	for(int j=js; j<je; j++){
		if(tpred[j] == time) continue;
		const double s = time - t0[j];
		const double s1 = 1.5 * s;
		const double s2 = 2.0 * s;
		for(int k=0; k<3; k++){
			x   [j][k] = x0   [j][k] + s *(x0dot[j][k] + s *(f   [j][k] + s*(fdot[j][k])));
			xdot[j][k] = x0dot[j][k] + s2*(f    [j][k] + s1*(fdot[j][k]));
		}
		tpred[j] = time;
	}
}


/*
 * SSE variant cxvpred_sse and supporting types (dbl6, mask2) removed.
 * The plain cxvpred above performs the same calculation without
 * x86-specific intrinsics.
 */

extern "C" void cxvpred_(
		const int *js,
		const int *je,
		const double *time,
		const double t0   [],
		const double x0   [][3],
		const double x0dot[][3],
		const double f    [][3],
		const double fdot [][3],
		_out_ double x    [][3],
		_out_ double xdot [][3],
		_out_ double tpred[])
{
	int cjs =*js;
	cjs--;
	int cje = *je;
	cje += cje%2;
        cxvpred(cjs, cje, *time, t0, x0, x0dot, f, fdot, x, xdot, tpred);
}
