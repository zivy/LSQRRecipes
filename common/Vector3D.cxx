#include "Vector3D.h"

namespace lsqrRecipes {

Vector3D crossProduct(const Vector3D &left, const Vector3D &right)
{
  Vector3D result;
  result[0] = left[1] * right[2] - left[2] * right[1];
  result[1] = left[2] * right[0] - left[0] * right[2];
  result[2] = left[0] * right[1] - left[1] * right[0];
  return result;
}

} //namespace lsqrRecipes
