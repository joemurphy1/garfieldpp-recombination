#ifndef G_GPUINTERFACE_H
#define G_GPUINTERFACE_H

// in theory, GPU code is faster with floats rather than doubles. This can be used to switch.
#define GPUFLOAT double

// Add this to do general optimisation of the GPU code. This will:
// * Remove 'unnecssary' checks
// * Disable GPU/CPU comparison code
// * Disable rounding code
//#define GPUOPTIMISE

//#define USEPRECALCRNG
#define MAXRNGDRAWS 14000

// Some constants for memory allocation
// should be dependent on hardware and dynamic, etc.
#define MAXCREATEDPARTICLES 2

#ifdef USEPRECALCRNG
#define MAXPARTICLES 250000
#define MAXSTACKSIZE 135000
#else
#define MAXPARTICLES 750000
#define MAXSTACKSIZE 500000
#endif

static const GPUFLOAT SmallGPU = 1.e-20;

#endif
