cmake_minimum_required(VERSION 2.8)

PROJECT( RANSAC_EXAMPLES )

option(USING_ITK_VNL "Build using ITK version of VNL." OFF)


#
# we use vnl from the vxl package to do all of our linear algebra and nonlinear
# optimization algorithms, either through ITK or directly through VXL
#
if(USING_ITK_VNL)
  #
  #remove the VXL_DIR from the cache, we are working with the ITK version
  #of vnl
  #
  unset(VXL_DIR CACHE)

  #
  # Using vnl through ITK, need to find ITK
  #
  find_package(ITK REQUIRED)
  if(ITK_FOUND)
    include(${ITK_USE_FILE})
    set(LINK_LIBRARIES
        LSQRRecipes
        itkvnl
        itkvnl_algo
        itkv3p_netlib
        itkvcl
        itkv3p_lsqr
    )
  endif()
else()
  #
  #remove the ITK_DIR from the cache, we are working with the VXL version of
  #vnl
  #
  unset(ITK_DIR CACHE)

  #
  # Using vnl through VXL, need to find VXL
  #
  find_package(VXL REQUIRED)
  if(VXL_FOUND)
    include(${VXL_CMAKE_DIR}/UseVXL.cmake)
    set(LINK_LIBRARIES
        LSQRRecipes
        vnl
        vnl_algo
        netlib
        v3p_netlib
        vcl
    )
  endif()
endif()


#
# get the settings required for using the LSQRRecipes library
#
find_package(LSQRRecipes)
if(LSQRRecipes_FOUND)
  include(${LSQRRecipes_USE_FILE})
endif()



#
# Example program showing the use of RANSAC for plane fitting
#
add_executable(planeEstimation planeEstimation.cxx)
target_link_libraries(planeEstimation ${LINK_LIBRARIES})


#
# Example program showing the use of RANSAC for sphere fitting
#
add_executable(sphereEstimation sphereEstimation.cxx)
target_link_libraries(sphereEstimation ${LINK_LIBRARIES})


#
# Example program showing the use of RANSAC for estimating the
# intersection of multiple 3D rays
#
add_executable(rayIntersectionEstimation rayIntersectionEstimation.cxx)
target_link_libraries(rayIntersectionEstimation ${LINK_LIBRARIES})


#
# Example program showing the use of RANSAC for estimating the
# calibration parameters of an US using a crosswire phantom
#
add_executable(crosswireUSCalibration crosswireUSCalibration.cxx)
target_link_libraries(crosswireUSCalibration ${LINK_LIBRARIES})


#
# Example program showing the use of RANSAC for estimating the
# calibration parameters of an US using a calibrated pointer
#
add_executable(pointerUSCalibration pointerUSCalibration.cxx)
target_link_libraries(pointerUSCalibration ${LINK_LIBRARIES})


#
# Example program showing the use of RANSAC for estimating the
# calibration parameters of an US using a plane phantom
#
add_executable(planeUSCalibration planeUSCalibration.cxx)
target_link_libraries(planeUSCalibration ${LINK_LIBRARIES})
