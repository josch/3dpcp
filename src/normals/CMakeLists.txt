IF(WITH_NORMALS)
  add_executable(calculateNormals calculate_normals.cc)

  target_link_libraries(calculateNormals scan ANN newmat ${Boost_LIBRARIES})
ENDIF(WITH_NORMALS)
