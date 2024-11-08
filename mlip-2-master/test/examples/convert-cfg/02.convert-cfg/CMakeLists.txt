set(TEST_NAME "02.convert-cfg")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/input.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/correctPOSCAR
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
  NAME ${TEST_NAME}
  COMMAND sh -c "($<TARGET_FILE:mlp> convert-cfg input.cfg ./out/POSCAR --output-format=vasp-poscar > /dev/null; diff correctPOSCAR ./out/POSCAR) > ./out/stdout.log 2>&1"
  WORKING_DIRECTORY ${TEST_DIR}
)