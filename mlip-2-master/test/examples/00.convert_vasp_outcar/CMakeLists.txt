set(TEST_NAME "convert_vasp_outcar")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/OUTCAR
  ${CMAKE_CURRENT_SOURCE_DIR}/sample_out
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
    NAME ${TEST_NAME}
    COMMAND sh -c "$<TARGET_FILE:mlp> convert-cfg --input-format=vasp-outcar OUTCAR ./out/relax.cfg > ./out/stdout.log"
    WORKING_DIRECTORY ${TEST_DIR}
)