IF(WITH_SEGMENTATION)
  add_executable(scan2segments scan2segments.cc)

  FIND_PACKAGE(OpenCV REQUIRED)

  target_link_libraries(scan2segments scan ANN fbr_cv_io fbr_panorama fbr_feature fbr_feature_matcher fbr_registration ${Boost_LIBRARIES} ${OpenCV_LIBS})


  set(CXX_FLAGS "-O2")
  set(FH_SEGMENTATION_SOURCES
	  main.cc
	  Options.cc
	  SRI.cpp
	  FHGraph.cpp
	  Timer.cpp

	  disjoint-set.cpp
	  segment-graph.cpp
	  segment-image.cpp
	  filter.cpp
	  convolve.cpp
  )

  set(FH_OUTPUT_SOURCES
	  fh_output.cpp
	  SRI.cpp

	  disjoint-set.cpp
	  segment-graph.cpp
	  segment-image.cpp
	  filter.cpp
	  convolve.cpp
  )

  add_executable(fh_segmentation ${FH_SEGMENTATION_SOURCES})

  target_link_libraries(fh_segmentation dl scanlib newmat ANN sparse scan_io_rxp scan_io_uosr scan_io_velodyne)
  
  add_executable(fh_output ${FH_OUTPUT_SOURCES})
  target_link_libraries(fh_output dl scanlib newmat sparse ANN scan_io_rxp scan_io_uosr scan_io_velodyne)

  enable_testing()
  add_executable(test_fhgraph FHGraph.cpp Timer.cpp test_fhgraph.cpp)
  target_link_libraries(test_fhgraph ANN dl)
   
  add_test(FHGraph ${CMAKE_BINARY_DIR}/bin/test_fhgraph)

ENDIF(WITH_SEGMENTATION)
