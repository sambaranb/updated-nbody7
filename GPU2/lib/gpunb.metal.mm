#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <cstdio>
#include <cmath>
#include <cassert>

#define NTHREAD 128

#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#else
static double get_wtime(){
	return 0.0;
}
#endif

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

// Metal state
static id<MTLDevice>               mtlDevice;
static id<MTLCommandQueue>         mtlQueue;
static id<MTLComputePipelineState> mtlPipeline;
static bool mtlInitialized = false;

// Persistent j-particle buffer
static id<MTLBuffer> jpBuf;
static Jparticle    *jpHost;
static int nbody, nbodymax;

// Reusable per-call buffers (avoid allocation every GPUNB_regf call)
static id<MTLBuffer> ipBuf;
static id<MTLBuffer> foBuf;
static id<MTLBuffer> nblBuf;
static int ipBufCount, nblBufNbmax;

static double time_send, time_grav;
static long long numInter;

static void metal_init(){
	if(mtlInitialized) return;

	mtlDevice = MTLCreateSystemDefaultDevice();
	if(!mtlDevice){
		fprintf(stderr, "gpunb: Metal is not supported on this device.\n");
		return;
	}
	mtlQueue = [mtlDevice newCommandQueue];

	NSError *error = nil;
	NSString *shaderFile = @"gpunb.metallib";
	NSArray *searchPaths = @[
		@".",
		[[NSString stringWithUTF8String:__FILE__] stringByDeletingLastPathComponent],
		@"./lib"
	];
	id<MTLLibrary> library = nil;
	for(NSString *dir in searchPaths){
		NSString *path = [dir stringByAppendingPathComponent:shaderFile];
		NSURL *url = [NSURL fileURLWithPath:path];
		library = [mtlDevice newLibraryWithURL:url error:&error];
		if(library) break;
	}
	if(!library){
		fprintf(stderr, "gpunb: failed to load %s: %s\n",
			[shaderFile UTF8String], [[error localizedDescription] UTF8String]);
		return;
	}
	id<MTLFunction> func = [library newFunctionWithName:@"gravity_kernel"];
	if(!func){
		fprintf(stderr, "gpunb: kernel 'gravity_kernel' not found.\n");
		return;
	}
	mtlPipeline = [mtlDevice newComputePipelineStateWithFunction:func error:&error];
	if(!mtlPipeline){
		fprintf(stderr, "gpunb: failed to create pipeline: %s\n",
			[[error localizedDescription] UTF8String]);
		return;
	}
	fprintf(stderr, "***********************\n");
	fprintf(stderr, "gpunb: Metal initialized on %s\n",
		[[mtlDevice name] UTF8String]);
	fprintf(stderr, "***********************\n");
	mtlInitialized = true;
}

void GPUNB_open(int nbmax){
	time_send = time_grav = 0.0;
	numInter = 0;
	nbodymax = nbmax;

	metal_init();
	if(!mtlInitialized) return;

	// Allocate persistent j-particle buffer (padded to NTHREAD)
	int ng = NTHREAD * ((nbmax + NTHREAD - 1) / NTHREAD);
	jpBuf = [mtlDevice newBufferWithLength:ng * sizeof(Jparticle)
		options:MTLResourceStorageModeShared];
	jpHost = (Jparticle *)[jpBuf contents];

	// Zero-fill the padding region
	for(int i=0; i<ng; i++){
		Jparticle p = {};
		jpHost[i] = p;
	}

#ifdef PROFILE
	fprintf(stderr, "***********************\n");
	fprintf(stderr, "Opened NBODY6/Metal library\n");
	fprintf(stderr, "nbmax = %d\n", nbmax);
	fprintf(stderr, "***********************\n");
#endif
}

void GPUNB_close(){
	jpBuf = nil;
	jpHost = nullptr;
	ipBuf = nil;
	foBuf = nil;
	nblBuf = nil;
	ipBufCount = 0;
	nblBufNbmax = 0;
	nbody = 0;
	nbodymax = 0;

#ifdef PROFILE
	fprintf(stderr, "***********************\n");
	fprintf(stderr, "Closed NBODY6/Metal library\n");
	fprintf(stderr, "time send   : %f sec\n", time_send);
	fprintf(stderr, "time grav   : %f sec\n", time_grav);
	if(time_grav > 0.0)
		fprintf(stderr, "%f Gflops (gravity part only)\n", 60.e-9 * numInter / time_grav);
	fprintf(stderr, "***********************\n");
#endif
}

void GPUNB_send(
		int nj,
		double mj[],
		double xj[][3],
		double vj[][3]){
	nbody = nj;
	assert(nbody <= nbodymax);
	time_send -= get_wtime();
	for(int j=0; j<nj; j++){
		jpHost[j].x[0] = (float)xj[j][0];
		jpHost[j].x[1] = (float)xj[j][1];
		jpHost[j].x[2] = (float)xj[j][2];
		jpHost[j].m    = (float)mj[j];
		jpHost[j].v[0] = (float)vj[j][0];
		jpHost[j].v[1] = (float)vj[j][1];
		jpHost[j].v[2] = (float)vj[j][2];
		jpHost[j].pad  = 0.f;
	}
	time_send += get_wtime();
}

void GPUNB_regf(
		int ni,
		double h2d[],
		double dtr[],
		double xid[][3],
		double vid[][3],
		double acc[][3],
		double jrk[][3],
		double pot[],
		int lmax,
		int nbmax,
		int *listbase){
	if(!mtlInitialized){
		fprintf(stderr, "gpunb: Metal not available.\n");
		return;
	}

	time_grav -= get_wtime();
	numInter += (long long)ni * nbody;

	int ng = NTHREAD * ((ni + NTHREAD - 1) / NTHREAD);
	int nnbmax = nbmax;

	// Check what needs reallocation before changing any state
	bool need_ip  = (!ipBuf  || ipBufCount < ng);
	bool need_nbl = (!nblBuf || nblBufNbmax < nnbmax || ipBufCount < ng);

	if(need_ip){
		ipBuf = [mtlDevice newBufferWithLength:ng * sizeof(Iparticle)
			options:MTLResourceStorageModeShared];
		foBuf = [mtlDevice newBufferWithLength:ng * sizeof(Force)
			options:MTLResourceStorageModeShared];
		ipBufCount = ng;
	}
	if(need_nbl){
		nblBuf = [mtlDevice newBufferWithLength:(long)ng * nnbmax * sizeof(int)
			options:MTLResourceStorageModeShared];
		nblBufNbmax = nnbmax;
	}

	Iparticle *ipHost = (Iparticle *)[ipBuf contents];
	for(int i=0; i<ni; i++){
		ipHost[i].x[0] = (float)xid[i][0];
		ipHost[i].x[1] = (float)xid[i][1];
		ipHost[i].x[2] = (float)xid[i][2];
		ipHost[i].h2   = (float)h2d[i];
		ipHost[i].v[0] = (float)vid[i][0];
		ipHost[i].v[1] = (float)vid[i][1];
		ipHost[i].v[2] = (float)vid[i][2];
		ipHost[i].dtr  = (float)dtr[i];
	}
	// Zero-fill padding
	Iparticle pzero = {};
	for(int i=ni; i<ng; i++){
		ipHost[i] = pzero;
	}

	// Dispatch inside autoreleasepool to free command buffer objects
	@autoreleasepool {
		id<MTLCommandBuffer> cmdBuf = [mtlQueue commandBuffer];
		id<MTLComputeCommandEncoder> encoder = [cmdBuf computeCommandEncoder];
		[encoder setComputePipelineState:mtlPipeline];
		[encoder setBuffer:jpBuf  offset:0 atIndex:0];
		[encoder setBuffer:ipBuf  offset:0 atIndex:1];
		[encoder setBuffer:foBuf  offset:0 atIndex:2];
		[encoder setBuffer:nblBuf offset:0 atIndex:3];
		[encoder setBytes:&nbody  length:sizeof(int) atIndex:4];
		[encoder setBytes:&nnbmax length:sizeof(int) atIndex:5];
		MTLSize gridSize  = MTLSizeMake(ng, 1, 1);
		MTLSize groupSize = MTLSizeMake(NTHREAD, 1, 1);
		[encoder dispatchThreads:gridSize threadsPerThreadgroup:groupSize];
		[encoder endEncoding];
		[cmdBuf commit];
		[cmdBuf waitUntilCompleted];
	}

	// Read results
	Force *foHost = (Force *)[foBuf contents];
	int   *nblHost = (int *)[nblBuf contents];

	for(int i=0; i<ni; i++){
		acc[i][0] = foHost[i].acc[0];
		acc[i][1] = foHost[i].acc[1];
		acc[i][2] = foHost[i].acc[2];
		jrk[i][0] = foHost[i].jrk[0];
		jrk[i][1] = foHost[i].jrk[1];
		jrk[i][2] = foHost[i].jrk[2];
		pot[i]    = foHost[i].pot;

		// Copy neighbor list in Fortran style
		int *nnbp    = listbase + lmax * i;
		int *nblistp = nnbp + 1;
		int nnb = foHost[i].nnb;
		if(nnb < 0 || nnb > nbmax){
			*nnbp = -1;
		}else{
			*nnbp = nnb;
			int off = i * nnbmax;
			for(int k=0; k<nnb; k++){
				nblistp[k] = nblHost[off + k];
			}
		}
	}

	time_grav += get_wtime();
}

extern "C" {
	void gpunb_devinit_(){
		metal_init();
	}
	void gpunb_open_(int *nbmax){
		GPUNB_open(*nbmax);
	}
	void gpunb_close_(){
		GPUNB_close();
	}
	void gpunb_send_(
			int *nj,
			double mj[],
			double xj[][3],
			double vj[][3]){
		GPUNB_send(*nj, mj, xj, vj);
	}
	void gpunb_regf_(
			int *ni,
			double h2[],
			double dtr[],
			double xi[][3],
			double vi[][3],
			double acc[][3],
			double jrk[][3],
			double pot[],
			int *lmax,
			int *nbmax,
			int *list){
		GPUNB_regf(*ni, h2, dtr, xi, vi, acc, jrk, pot, *lmax, *nbmax, list);
	}
}
