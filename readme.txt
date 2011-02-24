This repository contains implementations of algorithms for least squares
estimation of various parametric objects. In addition, a generic implementation of
the RANdom SAmple Consensus (RANSAC) is provided. This adds the ability to perform
robust estimation without integrating the associated logic into the least
squares logic.

The context within which these algorithms were used is image-guided interventions.

The parametric objects include:
1. (hyper) sphere.
2. (hyper) plane.
3. Intersection point of 3D rays.
4. US calibration matrix via (a) crosswire phantom (b) calibrated pointer
   (c) plane phantom.

Example and testing programs are provided.

Requirements:

Configuration is done via the Cmake tool (http://www.cmake.org).

The VNL subcomponent from the VXL (http://vxl.sourceforge.net/) toolkit is required
for linear algebra and nonlinear estimation methods. This set of libraries is also
available as part of ITK (http://www.itk.org/). The LSQRRecipes library can be
configured to use either of these distributions through a Cmake flag
(USING_ITK_VNL)
