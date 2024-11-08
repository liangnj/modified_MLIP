set(TEST_NAME "05.convert-cfg")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/input.bin.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_POSCAR0
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_POSCAR1
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_POSCAR2
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_POSCAR3
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
  NAME ${TEST_NAME}
  COMMAND sh -c "($<TARGET_FILE:mlp> convert-cfg input.bin.cfg ./out/POSCAR --input-format=bin --output-format=vasp-poscar > /dev/null; diff correct_POSCAR0 ./out/POSCAR0; diff correct_POSCAR1 ./out/POSCAR1; diff correct_POSCAR2 ./out/POSCAR2; diff correct_POSCAR3 ./out/POSCAR3) > ./out/stdout.log 2>&1"
  WORKING_DIRECTORY ${TEST_DIR}
)