This directory contains examples showing the usage of the RANSAC
algorithm to estimate various parametric objects.

The following examples also output 3D scenes:
1. lineEstimation
2. planeEstimation
3. rayIntersection
4. sphereEstimation
5. pivotCalibration

These 3D scenes visually illustrate the results of using a least
squares solution as compared to wrapping it with the RANSAC
algorithm. The output file format is OpenInventor. You can display
it using the viewer found in the coin3DViewer directory
(win32 command line application).

On linux (ubuntu) you can use the ivview to view open inventor
files, unfortunately it currently has a bug and does not run, but
it should be worth trying.

The following examples have data files associated with them:
1. linearEquationSystemSolver.
2. pivotCalibration.

Data/augmentedMatrixWithOutliers.txt
  This file contains an augmented matrix [A|b] with some rows being
  outliers. The example shows that we can robustly solve the
  overdetermined equation system Ax=b even with a significant amount
  of outliers.
Data/pivotCalibrationDataWithOutliers.txt
  This file contains a set of rigid transformation acquired while
  pivoting a tracked tool around a fixed point with about 30% outliers (arbitrary
  motion of the tool). Each line in the file is of the form x y z qx qy qz qs,
  where the rotation is represented by a unit quaternion.
