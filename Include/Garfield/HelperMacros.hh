#if defined(GARFIELD_CLASS_NAME)
  #undef GARFIELD_CLASS_NAME
#endif

#if defined(__GPUCOMPILE__)
  #define GARFIELD_CLASS_NAME(name) name ## GPU
#else 
  #define GARFIELD_CLASS_NAME(name) name
  //#define __device__
  //#define __host__
#endif
