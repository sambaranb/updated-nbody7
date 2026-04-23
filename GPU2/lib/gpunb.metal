#include <metal_stdlib>
using namespace metal;

#define NTHREAD 128

struct Jparticle{
	float x[3];
	float m;
	float v[3];
	float pad;
};

struct Iparticle{
	float x[3];
	float h2;
	float v[3];
	float dtr;
};

struct Force{
	float acc[3];
	float pot;
	float jrk[3];
	int   nnb;
};

kernel void gravity_kernel(
		device const Jparticle *jp  [[buffer(0)]],
		device const Iparticle *ip  [[buffer(1)]],
		device Force           *fo  [[buffer(2)]],
		device int             *nbl [[buffer(3)]],
		constant int &nbody         [[buffer(4)]],
		constant int &nnbmax        [[buffer(5)]],
		uint gid [[thread_position_in_grid]],
		uint lid [[thread_position_in_threadgroup]])
{
	threadgroup Jparticle jpbuf[NTHREAD];

	Iparticle myip = ip[gid];
	float ax = 0.f, ay = 0.f, az = 0.f;
	float jx = 0.f, jy = 0.f, jz = 0.f;
	float poti = 0.f;
	int nnb = 0;
	int nboff = gid * nnbmax;

	for(int j=0; j<nbody; j+=NTHREAD){
		threadgroup_barrier(mem_flags::mem_threadgroup);
		int jload = j + lid;
		if(jload < nbody)
			jpbuf[lid] = jp[jload];
		threadgroup_barrier(mem_flags::mem_threadgroup);

		int jmax = min(NTHREAD, nbody - j);
		for(int jj=0; jj<jmax; jj++){
			Jparticle jpart = jpbuf[jj];
			float dx  = jpart.x[0] - myip.x[0];
			float dy  = jpart.x[1] - myip.x[1];
			float dz  = jpart.x[2] - myip.x[2];
			float dvx = jpart.v[0] - myip.v[0];
			float dvy = jpart.v[1] - myip.v[1];
			float dvz = jpart.v[2] - myip.v[2];

			float dxp = dx + myip.dtr * dvx;
			float dyp = dy + myip.dtr * dvy;
			float dzp = dz + myip.dtr * dvz;

			float r2  = dx*dx   + dy*dy   + dz*dz;
			float r2p = dxp*dxp + dyp*dyp + dzp*dzp;
			float mh2 = jpart.m * myip.h2;

			float r2min = min(r2, r2p);
			bool is_nb = (r2min < mh2);

			float rinv1 = rsqrt(r2);
			if(is_nb){
				if(nnb < nnbmax)
					nbl[nboff + nnb] = j + jj;
				nnb++;
				rinv1 = 0.f;
			}

			float rv = dx*dvx + dy*dvy + dz*dvz;
			float rinv2 = rinv1 * rinv1;
			float mrinv1 = jpart.m * rinv1;
			float mrinv3 = mrinv1 * rinv2;
			rv *= -3.f * rinv2;

			poti += mrinv1;
			ax += mrinv3 * dx;
			ay += mrinv3 * dy;
			az += mrinv3 * dz;
			jx += mrinv3 * (dvx + rv * dx);
			jy += mrinv3 * (dvy + rv * dy);
			jz += mrinv3 * (dvz + rv * dz);
		}
	}

	Force result;
	result.acc[0] = ax;
	result.acc[1] = ay;
	result.acc[2] = az;
	result.pot    = poti;
	result.jrk[0] = jx;
	result.jrk[1] = jy;
	result.jrk[2] = jz;
	result.nnb    = (nnb <= nnbmax) ? nnb : -1;
	fo[gid] = result;
}
