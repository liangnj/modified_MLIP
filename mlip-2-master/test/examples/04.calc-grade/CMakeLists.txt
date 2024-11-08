set(TEST_NAME "calc-grade")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/pot.mtp
  ${CMAKE_CURRENT_SOURCE_DIR}/test.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/train.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/sample_out
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
    NAME ${TEST_NAME}
    COMMAND sh -c "$<TARGET_FILE:mlp> calc-grade pot.mtp train.cfg test.cfg ./out/graded.cfg --als-filename=./out/state.als > ./out/stdout.log"
    WORKING_DIRECTORY ${TEST_DIR}
)