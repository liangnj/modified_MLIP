set(TEST_NAME "14.convert-cfg")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/OUTCAR_relax
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_output.cfg
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
  NAME ${TEST_NAME}
  COMMAND sh -c "($<TARGET_FILE:mlp> convert-cfg OUTCAR_relax ./out/output.cfg --input-format=vasp-outcar --output-format=txt > /dev/null 2>&1; diff correct_output.cfg ./out/output.cfg) > ./out/stdout.log 2>&1"
  WORKING_DIRECTORY ${TEST_DIR}
)