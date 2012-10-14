IF(WITH_SEGMENTATION)
  add_executable(scan2segments scan2segments.cc)

  FIND_PACKAGE(OpenCV REQUIRED)

  target_link_libraries(scan2segments scan ANN fbr_cv_io fbr_panorama fbr_feature fbr_feature_matcher fbr_registration ${Boost_LIBRARIES} ${OpenCV_LIBS})


  set(CXX_FLAGS "-O2")
  set(FH_SEGMENTATION_SOURCES
	  fh_segmentation.cc
	  FHGraph.cpp

	  disjoint-set.cpp
	  segment-graph.cpp
	  segment-image.cpp
	  filter.cpp
	  convolve.cpp
  )

  set(FH_OUTPUT_SOURCES
	  fh_output.cpp

	  disjoint-set.cpp
	  segment-graph.cpp
	  segment-image.cpp
	  filter.cpp
	  convolve.cpp
  )


  set(SEG_LIBS dl scan newmat ANN sparse scanclient scanio)
  
  add_executable(fh_segmentation ${FH_SEGMENTATION_SOURCES})
  target_link_libraries(fh_segmentation ${SEG_LIBS}) 
  
  add_executable(fh_output ${FH_OUTPUT_SOURCES})
  target_link_libraries(fh_output ${SEG_LIBS}) 

  enable_testing()
  add_executable(test_fhgraph FHGraph.cpp test_fhgraph.cpp)
  target_link_libraries(test_fhgraph ${SEG_LIBS})
   
  add_test(FHGraph ${CMAKE_BINARY_DIR}/bin/test_fhgraph)

ENDIF(WITH_SEGMENTATION)
