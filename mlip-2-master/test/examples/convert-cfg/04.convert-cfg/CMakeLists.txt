set(TEST_NAME "04.convert-cfg")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/input.bin.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_output.cfg
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
  NAME ${TEST_NAME}
  COMMAND sh -c "($<TARGET_FILE:mlp> convert-cfg input.bin.cfg ./out/output.cfg --input-format=bin > /dev/null; diff correct_output.cfg ./out/output.cfg) > ./out/stdout.log 2>&1"
  WORKING_DIRECTORY ${TEST_DIR}
)