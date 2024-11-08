set(TEST_NAME "relax-active-learning")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/pot.mtp
  ${CMAKE_CURRENT_SOURCE_DIR}/mlip.ini
  ${CMAKE_CURRENT_SOURCE_DIR}/to-relax.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/train.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/state.als
  ${CMAKE_CURRENT_SOURCE_DIR}/sample_out
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
    NAME ${TEST_NAME}
    COMMAND sh -c "$<TARGET_FILE:mlp> relax mlip.ini --cfg-filename=to-relax.cfg --save-relaxed=./out/relaxed.cfg --save-unrelaxed=./out/unrelaxed.cfg > ./out/stdout.log"
    WORKING_DIRECTORY ${TEST_DIR}
)