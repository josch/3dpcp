IF(WITH_MIRRORDETECTION)
  add_executable(mirrordetection mirrordetection.cc)

  FIND_PACKAGE(OpenCV REQUIRED)

  target_link_libraries(mirrordetection scan ANN ${Boost_LIBRARIES} ${OpenCV_LIBS})
ENDIF(WITH_MIRRORDETECTION)
