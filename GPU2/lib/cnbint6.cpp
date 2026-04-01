#include <algorithm>
#include <cmath>
#include <cassert>

void cnbint(
		int i,
		const double pos[][3],
		const double vel[][3],
		const double mass[],
		int nnb,
		int list[],
		double f[3],
		double fdot[3]){
	const int NBMAX = 512;

	float xbuf[NBMAX];
	float ybuf[NBMAX];
	float zbuf[NBMAX];
	float vxbuf[NBMAX];
	float vybuf[NBMAX];
	float vzbuf[NBMAX];
	float mbuf[NBMAX];

	double xi = pos[i][0];
	double yi = pos[i][1];
	double zi = pos[i][2];
	float vxi = vel[i][0];
	float vyi = vel[i][1];
	float vzi = vel[i][2];
	double ax = 0.0, ay = 0.0, az = 0.0;
	float jx = 0.0f, jy = 0.0f, jz = 0.0f;
	for(int koff=0; koff<nnb; koff+=NBMAX){
		int nk = std::min(nnb-koff, NBMAX);
		for(int k=0; k<nk; k++){
			int j = list[k+koff];
			assert(j != i);
			double xj = pos[j][0];
			double yj = pos[j][1];
			double zj = pos[j][2];
			float vxj = vel[j][0];
			float vyj = vel[j][1];
			float vzj = vel[j][2];
			float mj = mass[j];
			xj -= xi;
			yj -= yi;
			zj -= zi;
			vxj -= vxi;
			vyj -= vyi;
			vzj -= vzi;
			xbuf[k] = xj;
			ybuf[k] = yj;
			zbuf[k] = zj;
			vxbuf[k] = vxj;
			vybuf[k] = vyj;
			vzbuf[k] = vzj;
			mbuf[k] = mj;
		}

		for(int k=0; k<nk; k++){
			float dx = xbuf[k];
			float dy = ybuf[k];
			float dz = zbuf[k];
			float dvx = vxbuf[k];
			float dvy = vybuf[k];
			float dvz = vzbuf[k];
			float mj = mbuf[k];

			float r2 = dx*dx + dy*dy + dz*dz;
			float rv = dx*dvx + dy*dvy + dz*dvz;
			rv *= -3.0f;
			float rinv1 = 1.0f / sqrtf(r2);
			float rinv2 = rinv1 * rinv1;
			rv *= rinv2;
			float rinv3 = mj * rinv1 * rinv2;

			dx *= rinv3; ax += (double)dx;
			dy *= rinv3; ay += (double)dy;
			dz *= rinv3; az += (double)dz;
			dvx *= rinv3; jx += dvx;
			dvy *= rinv3; jy += dvy;
			dvz *= rinv3; jz += dvz;
			dx *= rv; jx += dx;
			dy *= rv; jy += dy;
			dz *= rv; jz += dz;
		}
	}
	f[0] = ax;
	f[1] = ay;
	f[2] = az;
	fdot[0] = (double)jx;
	fdot[1] = (double)jy;
	fdot[2] = (double)jz;
	assert(f[0] == f[0]);
	assert(f[1] == f[1]);
	assert(f[2] == f[2]);
	assert(fdot[0] == fdot[0]);
	assert(fdot[1] == fdot[1]);
	assert(fdot[2] == fdot[2]);
}

extern "C" {
	void cnbint_(
		int *i,
		double pos[][3],
		double vel[][3],
		double mass[],
		int *nnb,
		int *nblist,
		double f[3],
		double fdot[3]){
		cnbint(*i, pos-1, vel-1, mass-1, *nnb-1, nblist, f, fdot);
	}
}
