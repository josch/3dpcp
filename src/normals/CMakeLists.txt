IF(WITH_NORMALS)
  FIND_PACKAGE(OpenCV REQUIRED)

  add_executable(calculateNormals calculate_normals.cc SRI.cc point.cc ../slam6d/fbr/fbr_global.cc ../slam6d/fbr/panorama.cc)

  target_link_libraries(calculateNormals scan ANN newmat fbr_feature fbr_cv_io fbr_panorama fbr_feature_matcher fbr_registration ${Boost_LIBRARIES} ${OpenCV_LIBS})
ENDIF(WITH_NORMALS)
