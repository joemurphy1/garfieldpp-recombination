#if defined(GARFIELD_CLASS_NAME)
  #undef GARFIELD_CLASS_NAME
#endif

#if defined(__GPUCOMPILE__)
  #define GARFIELD_CLASS_NAME(name) name ## GPU
  #undef __DEVICE__
  #define __DEVICE__ __device__
#else 
  #define GARFIELD_CLASS_NAME(name) name
  #undef __DEVICE__
  #define __DEVICE__
#endif
