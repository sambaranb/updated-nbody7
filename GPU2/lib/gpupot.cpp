#include <cmath>

struct float2{
	float x, y;
};
static inline float2 float2_split(double x){
	float2 ret;
	x *= (1<<16);
	double xi = (int)x;
	double xf = x - xi;
	ret.x = xi * (1./(1<<16));
	ret.y = xf * (1./(1<<16));
	return ret;
}

struct Particle{
	float2 pos[3];
	float mass;
	float pad;

	Particle(double x[3], double m){
		pos[0] = float2_split(x[0]);
		pos[1] = float2_split(x[1]);
		pos[2] = float2_split(x[2]);
		mass = (float)m;
	}
	Particle(){
		pos[0].x = pos[0].y = pos[1].x = pos[1].y = pos[2].x = pos[2].y = mass = pad = 0.f;
	}
};

void gpupot(
		int n,
		double m[],
		double x[][3],
		double pot[]){
	Particle *ptcl = new Particle[n];
	for(int i=0; i<n; i++){
		ptcl[i] = Particle(x[i], m[i]);
	}

#pragma omp parallel for
	for(int i=0; i<n; i++){
		float potH = 0.f;
		float potL = 0.f;
		float xiH = ptcl[i].pos[0].x;
		float yiH = ptcl[i].pos[1].x;
		float ziH = ptcl[i].pos[2].x;
		float xiL = ptcl[i].pos[0].y;
		float yiL = ptcl[i].pos[1].y;
		float ziL = ptcl[i].pos[2].y;
		for(int j=0; j<n; j++){
			float xjH = ptcl[j].pos[0].x;
			float xjL = ptcl[j].pos[0].y;
			float yjH = ptcl[j].pos[1].x;
			float yjL = ptcl[j].pos[1].y;
			float zjH = ptcl[j].pos[2].x;
			float zjL = ptcl[j].pos[2].y;
			float mj  = ptcl[j].mass;

			float dx = (xjH - xiH) + (xjL - xiL);
			float dy = (yjH - yiH) + (yjL - yiL);
			float dz = (zjH - ziH) + (zjL - ziL);
			float r2 = dx*dx + dy*dy + dz*dz;
			float rinv = (r2 > 0.f) ? 1.f / sqrtf(r2) : 0.f;
			rinv *= mj;

			float tmp = potH;
			potH += rinv;
			potL -= (potH - tmp) - rinv;
		}
		pot[i] = (double)potH + (double)potL;
	}

	delete [] ptcl;
}

/* Path 2 C1 (internal MPI): evaluate the potential only for the i-slice
   [i0 .. i0+ni-1] (0-based) against the FULL j-set, so each MPI rank can
   compute its contiguous share and Allgather the result. The j-loop is a
   verbatim copy of gpupot()'s inner loop -- same order, same compensated
   float summation -- so pot[i] is bit-identical to the full evaluation
   (the np1-vs-npK equivalence harness verifies this on every change).
   gpupot() itself is deliberately untouched: the serial path keeps the
   exact shipped instruction stream. */
void gpupot_range(
		int i0,
		int ni,
		int n,
		double m[],
		double x[][3],
		double pot[]){
	Particle *ptcl = new Particle[n];
	for(int i=0; i<n; i++){
		ptcl[i] = Particle(x[i], m[i]);
	}

#pragma omp parallel for
	for(int ii=0; ii<ni; ii++){
		const int i = i0 + ii;
		float potH = 0.f;
		float potL = 0.f;
		float xiH = ptcl[i].pos[0].x;
		float yiH = ptcl[i].pos[1].x;
		float ziH = ptcl[i].pos[2].x;
		float xiL = ptcl[i].pos[0].y;
		float yiL = ptcl[i].pos[1].y;
		float ziL = ptcl[i].pos[2].y;
		for(int j=0; j<n; j++){
			float xjH = ptcl[j].pos[0].x;
			float xjL = ptcl[j].pos[0].y;
			float yjH = ptcl[j].pos[1].x;
			float yjL = ptcl[j].pos[1].y;
			float zjH = ptcl[j].pos[2].x;
			float zjL = ptcl[j].pos[2].y;
			float mj  = ptcl[j].mass;

			float dx = (xjH - xiH) + (xjL - xiL);
			float dy = (yjH - yiH) + (yjL - yiL);
			float dz = (zjH - ziH) + (zjL - ziL);
			float r2 = dx*dx + dy*dy + dz*dz;
			float rinv = (r2 > 0.f) ? 1.f / sqrtf(r2) : 0.f;
			rinv *= mj;

			float tmp = potH;
			potH += rinv;
			potL -= (potH - tmp) - rinv;
		}
		pot[i] = (double)potH + (double)potL;
	}

	delete [] ptcl;
}

extern "C"{
	void gpupot_(
			int *n,
			double m[],
			double x[][3],
			double pot[]){
		gpupot(*n, m, x, pot);
	}
	/* Fortran binding; i0 is the 1-based first slice member. */
	void gpupot_range_(
			int *i0,
			int *ni,
			int *n,
			double m[],
			double x[][3],
			double pot[]){
		gpupot_range(*i0 - 1, *ni, *n, m, x, pot);
	}
}
