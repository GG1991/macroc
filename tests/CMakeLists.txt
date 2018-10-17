
set(CMAKE_BUILD_TYPE Debug)

# The test files I care only (create a list)
set(testsources
  test_dm_1.c
  test_dm_2.c
  )

# Iterate over the list above
foreach (testfile ${testsources})

  # Delete File extensions (test_i.cpp -> test_i)
  get_filename_component(testname ${testfile} NAME_WE)

  # Create the executable test_i and link.
  add_executable(${testname} ${testfile})
  target_link_libraries(${testname} ${PETSC_LIB})
endforeach ()
