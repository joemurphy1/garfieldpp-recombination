// Include this header if we're compiling with the GPU or this is the first time without
#if defined(__GPUCOMPILE__) || !defined(TETRAHEDRAL_TREE_H)

#if !defined(__GPUCOMPILE__) && !defined(TETRAHEDRAL_TREE_H)
#define TETRAHEDRAL_TREE_H
#endif

#include "Garfield/HelperMacros.hh"
#include "Garfield/Vector.hh"

#if defined(__GPUCOMPILE__)
  #include"GPUInterface.hh"
#else
  class TetrahedralTreeGPU;
#endif

#include <cstddef>
#include <vector>
#include <utility>

namespace Garfield
{

#ifdef __GPUCOMPILE__
  using Vec3 =  Vec3Impl<cuda_t>;
#else
  using Vec3 =  Vec3Impl<double>;
#endif

/**

\brief Helper class for searches in field maps.

This class stores the mesh nodes and elements in an Octree data
structure to optimize the element search operations

Author: Ali Sheharyar

Organization: Texas A&M University at Qatar

*/
class GARFIELD_CLASS_NAME(TetrahedralTree) {
 public:
 #ifdef __GPUCOMPILE__
  // Constructor
  GARFIELD_CLASS_NAME(TetrahedralTree)() = default;

  // Destructor
  ~GARFIELD_CLASS_NAME(TetrahedralTree)() = default;
 #else
  // Constructor
  GARFIELD_CLASS_NAME(TetrahedralTree)(const Vec3& origin, const Vec3& halfDimension);

  /// Destructor
  ~GARFIELD_CLASS_NAME(TetrahedralTree)();
#endif

  #ifndef __GPUCOMPILE__
  // Insert a mesh node (a vertex/point) to the tree
  void InsertMeshNode(Vec3 point, const int index);

  /// Insert a mesh element with given bounding box and index to the tree.
  void InsertMeshElement(const double bb[6], const int index);

  /// Create and initialise GPU Transfer class
  double CreateGPUTransferObject(TetrahedralTreeGPU *&tree_gpu);
  #endif

 private:
  static std::vector<int> emptyBlock;

  // Physical centre of this tree node.
  Vec3 m_origin;
  // Half the width/height/depth of this tree node. 
  Vec3 m_halfDimension;
  // Storing min and max points for convenience
  Vec3 m_min, m_max;  

  // The tree has up to eight children and can additionally store
  // a list of mesh nodes and mesh elements.
  // Pointers to child octants.
  GARFIELD_CLASS_NAME(TetrahedralTree)* children[8];  

  // Children follow a predictable pattern to make accesses simple.
  // Here, - means less than 'origin' in that dimension, + means greater than.
  // child:	0 1 2 3 4 5 6 7
  // x:     - - - - + + + +
  // y:     - - + + - - + +
  // z:     - + - + - + - +

  #ifndef __GPUCOMPILE__
  std::vector<std::pair<Vec3, int> > nodes;
  #endif

  #ifdef __GPUCOMPILE__
  int* elements{nullptr};
  int numelements{0};
  #else
  std::vector<int> elements;
  #endif

  static const size_t BlockCapacity = 10;

  #ifndef __GPUCOMPILE__
  // Check if the given box overlaps with this tree node.
  bool DoesBoxOverlap(const double bb[6]) const;
  #endif
  // Check if this tree node is a leaf or intermediate node.
  __DEVICE__ bool IsLeafNode() const;


public:
  // Get all tetrahedra linked to a block corresponding to the given point
  #ifdef __GPUCOMPILE__
  __device__ void GetElementsInBlock(const Vec3& point, const int *&tet_list_elems, int &num_elems) const;
  #else
  const std::vector<int>& GetElementsInBlock(const Vec3& point) const;
  #endif

private:
__DEVICE__ int GetOctantContainingPoint(const Vec3& point) const;
  // Get a block containing the input point
  __DEVICE__ const GARFIELD_CLASS_NAME(TetrahedralTree)* GetBlockFromPoint(const Vec3& point) const;
  // A helper function used by the function above.
  // Called recursively on the child nodes.
  __DEVICE__ const GARFIELD_CLASS_NAME(TetrahedralTree)* GetBlockFromPointHelper(const Vec3& point) const;  

#ifdef __GPUCOMPILE__
  friend class TetrahedralTree;
#endif
};

}

#endif
