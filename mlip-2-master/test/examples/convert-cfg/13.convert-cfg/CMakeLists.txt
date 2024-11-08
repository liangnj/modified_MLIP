set(TEST_NAME "13.convert-cfg")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/OUTCAR
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_output.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/output.cfg
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})
file(COPY ${TEST_DIR}/output.cfg DESTINATION ${TEST_DIR}/out)

add_test(
  NAME ${TEST_NAME}
  COMMAND sh -c "($<TARGET_FILE:mlp> convert-cfg OUTCAR ./out/output.cfg --input-format=vasp-outcar --append > /dev/null; diff correct_output.cfg ./out/output.cfg) > ./out/stdout.log 2>&1"
  WORKING_DIRECTORY ${TEST_DIR}
)