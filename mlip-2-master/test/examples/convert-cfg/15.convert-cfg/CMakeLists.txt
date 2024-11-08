set(TEST_NAME "15.convert-cfg")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/sample_input.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_lammps.inp
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
  NAME ${TEST_NAME}
  COMMAND sh -c "($<TARGET_FILE:mlp> convert-cfg sample_input.cfg ./out/lammps.inp --output-format=lammps-datafile --input-format=txt > /dev/null 2>&1; diff correct_lammps.inp ./out/lammps.inp) > ./out/stdout.log 2>&1"
  WORKING_DIRECTORY ${TEST_DIR}
)