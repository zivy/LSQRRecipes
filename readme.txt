This repository contains implementations of algorithms for least squares
estimation of various parametric objects. In addition, a generic implementation of
the RANdom SAmple Consensus (RANSAC) is provided. This adds the ability to perform
robust estimation without integrating the associated logic into the least
squares logic.

The context within which these algorithms were used is image-guided interventions.

Example and testing programs are provided.

The parametric objects include:
1. hypersphere (2D - circle, 3D - sphere, kD - sphere).
2. hyperplane (2D - line, 3D - plane, kD - plane).
3. intersection point of 3D rays.
4. US calibration matrix via (a) crosswire phantom (b) calibrated pointer
   (c) plane phantom.
5. absolute orientation, rigid transformation, relating corresponding point pairs
   in two coordinate systems [unit quaternion solution due to B.K.P. Horn].
6. kD line (line in any dimension).
7. 2D line - doesn't require the VNL library and can be easily incorporated
   into an existing project using a minimal subset of the LSQRRecipes library
   (see details below).
8. solution to a set of linear equations.
9. translations obtained via pivot calibration.
-------------------------------------------------------------------------------

Requirements:

Configuration is done via the Cmake tool (http://www.cmake.org).

The VNL subcomponent from the VXL (http://vxl.sourceforge.net/) toolkit is required
for linear algebra and nonlinear estimation methods. This set of libraries is also
available as part of ITK (http://www.itk.org/). The LSQRRecipes library can be
configured to use either of these distributions through a Cmake flag
(USING_ITK_VNL)

-------------------------------------------------------------------------------

2D line
--------
If all you need is a robust estimation of a 2D line and you don't want to use a
linear algebra package for computing eigenvectors you can directly copy the
following files into your project:
0. copyright.h
1. ParametersEstimator.h
2. Line2DParametersEstimator.{h,cxx}
3. RANSAC.{h,hxx}

You will need a 2D point class which has the square bracket operator
to access the point elements (Point2D.h). Either you already have such a class
or you can use the following minimal implementation, just copy paste:

#ifndef _POINT2D_H_
#define _POINT2D_H_

#include "copyright.h"

namespace lsqrRecipes {

class Point2D
{
public:
  Point2D(double x=0.0, double y=0.0) {this->data[0] = x; this->data[1] = y;}
  const double & operator[](int index) const {return this->data[index];}

private:
  double data[2];
};

} //namespace lsqrRecipes

#endif //_POINT2D_H_
-------------------------------------------------------------------------------
