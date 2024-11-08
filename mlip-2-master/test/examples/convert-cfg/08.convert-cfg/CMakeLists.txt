set(TEST_NAME "08.convert-cfg")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/OUTCAR
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_output.bin.cfg
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
  NAME ${TEST_NAME}
  COMMAND sh -c "($<TARGET_FILE:mlp> convert-cfg OUTCAR ./out/output.bin.cfg --input-format=vasp-outcar --output-format=bin > /dev/null; diff correct_output.bin.cfg ./out/output.bin.cfg) > ./out/stdout.log 2>&1"
  WORKING_DIRECTORY ${TEST_DIR}
)