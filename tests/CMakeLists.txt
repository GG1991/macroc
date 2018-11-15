
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

add_test(NAME test_coupling_small_1 COMMAND macroc -da_grid_x 3 -da_grid_y 3 -da_grid_z 3 -ts 5 > test_coupling_small_1.out)
add_test(NAME test_coupling_small_2 COMMAND macroc -da_grid_x 4 -da_grid_y 4 -da_grid_z 4 -ts 5 > test_coupling_small_2.out)
add_test(NAME test_coupling_small_3 COMMAND macroc -da_grid_x 5 -da_grid_y 2 -da_grid_z 2 -ts 5 > test_coupling_small_3.out)
