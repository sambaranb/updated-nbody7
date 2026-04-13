#include <metal_stdlib>
using namespace metal;

#define NTHREAD 128

struct Particle{
	float2 pos[3];
	float mass;
	float pad;
};

static inline float2 float2_accum(float2 acc, float x){
	float tmp = acc.x + x;
	acc.y -= (tmp - acc.x) - x;
	acc.x = tmp;
	return acc;
}

static inline float2 float2_regularize(float2 acc){
	float tmp = acc.x + acc.y;
	acc.y = acc.y - (tmp - acc.x);
	acc.x = tmp;
	return acc;
}

kernel void pot_kernel(
		device const Particle *ptcl [[buffer(0)]],
		device float2 *phi          [[buffer(1)]],
		constant int &n             [[buffer(2)]],
		uint gid [[thread_position_in_grid]],
		uint lid [[thread_position_in_threadgroup]])
{
	threadgroup Particle jpbuf[NTHREAD];
	Particle ip = ptcl[gid];
	float2 phii = float2(0.f, 0.f);
	for(int j=0; j<n; j+=NTHREAD){
		threadgroup_barrier(mem_flags::mem_threadgroup);
		jpbuf[lid] = ptcl[j + lid];
		threadgroup_barrier(mem_flags::mem_threadgroup);
		for(int jj=0; jj<NTHREAD; jj++){
			Particle jp = jpbuf[jj];
			float dx = (jp.pos[0].x - ip.pos[0].x) + (jp.pos[0].y - ip.pos[0].y);
			float dy = (jp.pos[1].x - ip.pos[1].x) + (jp.pos[1].y - ip.pos[1].y);
			float dz = (jp.pos[2].x - ip.pos[2].x) + (jp.pos[2].y - ip.pos[2].y);
			float r2 = dx*dx + dy*dy + dz*dz;
			float pij = jp.mass * rsqrt(r2);
			if(r2 > 0.f) phii = float2_accum(phii, pij);
		}
		phii = float2_regularize(phii);
	}
	phi[gid] = phii;
}
