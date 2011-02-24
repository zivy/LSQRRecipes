
add_executable(planeEstimationTest PlaneParametersEstimatorTest.cxx)
target_link_libraries(planeEstimationTest ${LINK_LIBRARIES})

add_executable(sphereEstimationTest SphereParametersEstimatorTest.cxx)
target_link_libraries(sphereEstimationTest ${LINK_LIBRARIES})

add_executable(rayIntersectionEstimationTest RayIntersectionParametersTest.cxx)
target_link_libraries(rayIntersectionEstimationTest ${LINK_LIBRARIES})

add_executable(singlePointUSCalibrationEstimationTest SinglePointTargetUSCalibrationParametersEstimatorTest.cxx)
target_link_libraries(singlePointUSCalibrationEstimationTest ${LINK_LIBRARIES})

add_executable(planePhantomUSCalibrationEstimationTest PlanePhantomUSCalibrationParametersEstimatorTest.cxx)
target_link_libraries(planePhantomUSCalibrationEstimationTest ${LINK_LIBRARIES})


add_test(PlaneEstimationTest planeEstimationTest)
add_test(SphereEstimationTest sphereEstimationTest)
add_test(RayIntersectionEstimationTest rayIntersectionEstimationTest)
add_test(SinglePointTargetUSCalibrationSimulatedTest singlePointUSCalibrationEstimationTest)
add_test(PlanePhantomUSCalibrationSimulatedTest planePhantomUSCalibrationEstimationTest)

if(LSQR_RECIPES_TESTING_DATA_ROOT)
  add_test(SinglePointTargetUSCalibrationPhantomExperimentalDataTest singlePointUSCalibrationEstimationTest
           ${LSQR_RECIPES_TESTING_DATA_ROOT}/crossWirePhantomTransformations.txt
           ${LSQR_RECIPES_TESTING_DATA_ROOT}/crossWirePhantom2DPoints.txt)
else()
    message(FATAL_ERROR "LSQR_RECIPES_TESTING_DATA_ROOT should be set to the testing data root directory.")
endif()
