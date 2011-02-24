#
# This is the configuration file for the LSQRRecipes (least-squares recipes)
# library. The library requires the vnl components from the vxl project
# (http://vxl.sourceforge.net/). The configuration supports either the use of vnl
# from a vxl installation or ITK's (www.itk.org) version of vnl, if installed.
# The choice of which version to use is indicated by the Cmake variable USING_ITK_VNL.
#

PROJECT( LSQRRecipes )

cmake_minimum_required(VERSION 2.8)

#
# LSQRRecipes version number.
#
set(LSQRRecipes_VERSION_MAJOR "0")
set(LSQRRecipes_VERSION_MINOR "1")
set(LSQRRecipes_VERSION_PATCH "0")

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
# include directories
#
set(LSQRRecipes_INCLUDE_DIRS
    ${LSQRRecipes_SOURCE_DIR}/copyright
    ${LSQRRecipes_SOURCE_DIR}/common
    ${LSQRRecipes_SOURCE_DIR}/parametersEstimators
)

include_directories(
  ${LSQRRecipes_INCLUDE_DIRS}
  )


set(LSQRRecipes_HDRS
   copyright/copyright.h
   common/Epsilon.h
   common/RandomNumberGenerator.h
   common/Point.h
   common/Point2D.h
   common/IPoint2D.h
   common/Point3D.h
   common/Vector.h
   common/Vector3D.h
   common/Ray3D.h
   common/Frame.h
   parametersEstimators/ParametersEstimator.h
   parametersEstimators/RANSAC.h
   parametersEstimators/RANSAC.txx
   parametersEstimators/PlaneParametersEstimator.h
   parametersEstimators/PlaneParametersEstimator.txx
   parametersEstimators/SphereParametersEstimator.h
   parametersEstimators/SphereParametersEstimator.txx
   parametersEstimators/RayIntersectionParametersEstimator.h
   parametersEstimators/SinglePointTargetUSCalibrationParametersEstimator.h
   parametersEstimators/PlanePhantomUSCalibrationParametersEstimator.h
)

set(LSQRRecipes_SRCS
   common/Vector3D.cxx
   common/Ray3D.cxx
   common/Frame.cxx
   parametersEstimators/RayIntersectionParametersEstimator.cxx
   parametersEstimators/SinglePointTargetUSCalibrationParametersEstimator.cxx
   parametersEstimators/PlanePhantomUSCalibrationParametersEstimator.cxx
)

#
# Actually create the LSQRRecipes library
#
add_library(LSQRRecipes ${LSQRRecipes_SRCS} ${LSQRRecipes_HDRS})
set_target_properties(LSQRRecipes PROPERTIES
                      VERSION ${LSQRRecipes_VERSION_MAJOR}.${LSQRRecipes_VERSION_MINOR}.${LSQRRecipes_VERSION_PATCH}
                      SOVERSION ${LSQRRecipes_VERSION_MAJOR}
)

################################################################################
#
# optional testing of the estimation classes
#
option(BUILD_TESTING "Build the Testing directory." ON)
if(BUILD_TESTING)
  enable_testing()
  set(LSQR_RECIPES_TESTING_DATA_ROOT ${LSQRRecipes_SOURCE_DIR}/Testing/Data)
  add_subdirectory(testing)
endif()

################################################################################
#
#settings for installation.
#

#used by the GenerateLSQRRecipesConfig.cmake file
set(PRE_INSTALL_DIR ${LSQRRecipes_BINARY_DIR}/preInstall)

#
# where to install the library
#
set(LSQRRecipes_INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "Directory where you want to install the LSQRRecipes library.")

#
# install the static library
#
install(TARGETS LSQRRecipes
        ARCHIVE DESTINATION ${LSQRRecipes_INSTALL_LIB_DIR}
        LIBRARY DESTINATION ${LSQRRecipes_INSTALL_LIB_DIR}
)

#
# where to install the headers
#
set(LSQRRecipes_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include CACHE PATH "Directory where you want to install the LSQRRecipes library headers.")

#
# install all of the headers
#
install(FILES ${LSQRRecipes_HDRS}
        DESTINATION ${LSQRRecipes_INSTALL_INCLUDE_DIR}
)

#
# copy USeLSQRRecipes.cmake to the library dir so that external
# projects have the correct settings
#
install(FILES ${LSQRRecipes_BINARY_DIR}/UseLSQRRecipes.cmake
              ${PRE_INSTALL_DIR}/LSQRRecipesConfig.cmake
        DESTINATION ${LSQRRecipes_INSTALL_LIB_DIR}
)

################################################################################
#
# settings so external projects can easily use LSQRRecipes
#

#
# create the LSQRRecipes.cmake file, using the GenerateLSQRRecipes configuration
#
include(${LSQRRecipes_SOURCE_DIR}/GenerateLSQRRecipesConfig.cmake)

#
#copy the UseLSQRRecipes.cmake file to the binary dir
#
configure_file(${LSQRRecipes_SOURCE_DIR}/UseLSQRRecipes.cmake.in
               ${LSQRRecipes_BINARY_DIR}/UseLSQRRecipes.cmake COPYONLY)
