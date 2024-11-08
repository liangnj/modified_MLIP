set(TEST_NAME "06.convert-cfg")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/input.bin.cfg
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
  NAME ${TEST_NAME}
  COMMAND sh -c "($<TARGET_FILE:mlp> convert-cfg input.bin.cfg ./out/output.bin.cfg --input-format=bin --output-format=bin > /dev/null; diff input.bin.cfg ./out/output.bin.cfg) > ./out/stdout.log 2>&1"
  WORKING_DIRECTORY ${TEST_DIR}
)