IF(WITH_SEGMENTATION)
  add_executable(scan2segments scan2segments.cc ../slam6d/fbr/fbr_global.cc)

  FIND_PACKAGE(OpenCV REQUIRED)

  target_link_libraries(scan2segments scan ANN fbr_feature fbr_cv_io fbr_panorama fbr_feature_matcher fbr_registration ${Boost_LIBRARIES} ${OpenCV_LIBS})
ENDIF(WITH_SEGMENTATION)
