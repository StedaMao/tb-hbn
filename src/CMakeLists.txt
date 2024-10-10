set(EXECUTABLE tb_hbn)

add_executable(${EXECUTABLE}
  main.f90
  )


set_target_properties(${EXECUTABLE}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )

