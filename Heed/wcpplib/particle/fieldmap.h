#ifndef FIELDMAP_H
#define FIELDMAP_H

#include "wcpplib/clhep_units/WPhysicalConstants.h"
#include "wcpplib/geometry/vec.h"

namespace Heed {

/// Retrieve electric and magnetic field.

class fieldmap {
 public:
  fieldmap() = default;

  virtual void evaluate(const point& /*pt*/, vec& efield, vec& bfield,
                vfloat& mrange) const {
    efield.x = bfield.x = 0.;
    efield.y = bfield.y = 0.;
    efield.z = bfield.z = 0.;
    mrange = DBL_MAX;
  }
  virtual bool inside(const point& /*pt*/) { return true; }

};
}

#endif
