#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <cstdio>
#include <cmath>

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

struct Particle{
	float pos[6]; // pos[0].x, pos[0].y, pos[1].x, pos[1].y, pos[2].x, pos[2].y
	float mass;
	float pad;
};

static Particle make_particle(double x[3], double m){
	Particle p;
	const int shift = 20;
	for(int k=0; k<3; k++){
		double xk = x[k] * (1<<shift);
		double xi = (int)xk;
		double xf = xk - xi;
		p.pos[2*k]   = xi * (1./(1<<shift));
		p.pos[2*k+1] = xf * (1./(1<<shift));
	}
	p.mass = (float)m;
	p.pad  = 0.f;
	return p;
}

static Particle make_particle_zero(){
	Particle p;
	for(int k=0; k<6; k++) p.pos[k] = 0.f;
	p.mass = 0.f;
	p.pad  = 0.f;
	return p;
}

static id<MTLDevice>              mtlDevice;
static id<MTLCommandQueue>        mtlQueue;
static id<MTLComputePipelineState> mtlPipeline;
static bool mtlInitialized = false;

// Reusable buffers (avoid allocation every gpupot call)
static id<MTLBuffer> ptclBuf;
static id<MTLBuffer> phiBuf;
static int ptclBufCount;

static void metal_init(){
	if(mtlInitialized) return;

	mtlDevice = MTLCreateSystemDefaultDevice();
	if(!mtlDevice){
		fprintf(stderr, "gpupot: Metal is not supported on this device.\n");
		return;
	}
	mtlQueue = [mtlDevice newCommandQueue];

	// Load pre-compiled metallib from same directory as execution
	NSError *error = nil;
	NSString *libPath = [[NSString stringWithUTF8String:__FILE__]
		stringByDeletingLastPathComponent];
	// Look for metallib next to the executable first, then in lib/
	NSString *metallib = @"gpupot.metallib";
	NSArray *searchPaths = @[
		@".",
		libPath,
		@"./lib"
	];
	id<MTLLibrary> library = nil;
	for(NSString *dir in searchPaths){
		NSString *path = [dir stringByAppendingPathComponent:metallib];
		NSURL *url = [NSURL fileURLWithPath:path];
		library = [mtlDevice newLibraryWithURL:url error:&error];
		if(library) break;
	}
	if(!library){
		fprintf(stderr, "gpupot: failed to load %s: %s\n",
			[metallib UTF8String], [[error localizedDescription] UTF8String]);
		return;
	}
	id<MTLFunction> func = [library newFunctionWithName:@"pot_kernel"];
	if(!func){
		fprintf(stderr, "gpupot: kernel function 'pot_kernel' not found.\n");
		return;
	}
	mtlPipeline = [mtlDevice newComputePipelineStateWithFunction:func error:&error];
	if(!mtlPipeline){
		fprintf(stderr, "gpupot: failed to create pipeline: %s\n",
			[[error localizedDescription] UTF8String]);
		return;
	}
	fprintf(stderr, "gpupot: Metal initialized on %s\n",
		[[mtlDevice name] UTF8String]);
	mtlInitialized = true;
}

void gpupot(
		int n,
		double m[],
		double x[][3],
		double pot[]){
	metal_init();
	if(!mtlInitialized){
		fprintf(stderr, "gpupot: Metal not available, cannot compute.\n");
		return;
	}

	double t0 = get_wtime();

	int ng = NTHREAD * ((n + NTHREAD - 1) / NTHREAD);

	// Reuse buffers if large enough, otherwise reallocate
	if(!ptclBuf || ptclBufCount < ng){
		ptclBuf = [mtlDevice newBufferWithLength:ng * sizeof(Particle)
			options:MTLResourceStorageModeShared];
		phiBuf = [mtlDevice newBufferWithLength:ng * sizeof(float) * 2
			options:MTLResourceStorageModeShared];
		ptclBufCount = ng;
	}

	Particle *ptcl = (Particle *)[ptclBuf contents];
	for(int i=0; i<n; i++){
		ptcl[i] = make_particle(x[i], m[i]);
	}
	for(int i=n; i<ng; i++){
		ptcl[i] = make_particle_zero();
	}

	// Encode and dispatch inside autoreleasepool
	@autoreleasepool {
		id<MTLCommandBuffer> cmdBuf = [mtlQueue commandBuffer];
		id<MTLComputeCommandEncoder> encoder = [cmdBuf computeCommandEncoder];
		[encoder setComputePipelineState:mtlPipeline];
		[encoder setBuffer:ptclBuf offset:0 atIndex:0];
		[encoder setBuffer:phiBuf  offset:0 atIndex:1];
		[encoder setBytes:&n length:sizeof(int) atIndex:2];
		MTLSize gridSize = MTLSizeMake(ng, 1, 1);
		MTLSize groupSize = MTLSizeMake(NTHREAD, 1, 1);
		[encoder dispatchThreads:gridSize threadsPerThreadgroup:groupSize];
		[encoder endEncoding];
		[cmdBuf commit];
		[cmdBuf waitUntilCompleted];
	}

	// Read results
	float *phi = (float *)[phiBuf contents];
	for(int i=0; i<n; i++){
		pot[i] = (double)phi[2*i] + (double)phi[2*i+1];
	}

	double t1 = get_wtime();
#ifdef PROFILE
	fprintf(stderr, "gpupot: %f sec (Metal, n=%d)\n", t1 - t0, n);
#endif
}

extern "C"{
	void gpupot_(
			int *n,
			double m[],
			double x[][3],
			double pot[]){
		gpupot(*n, m, x, pot);
	}
}
