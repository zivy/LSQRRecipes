#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include "copyright.h"
#include "Vector.h"

namespace lsqrRecipes {

typedef Vector<double,3> Vector3D;

Vector3D crossProduct(const Vector3D &left, const Vector3D &right);

} //namespace lsqrRecipes

#endif //_VECTOR3D_H_
