set(TEST_NAME "train")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/06.mtp
  ${CMAKE_CURRENT_SOURCE_DIR}/train.cfg
  ${CMAKE_CURRENT_SOURCE_DIR}/sample_out
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
    NAME ${TEST_NAME}
    COMMAND sh -c "$<TARGET_FILE:mlp> train 06.mtp train.cfg --trained-pot-name=./out/pot.mtp --max-iter=100 > ./out/stdout.log"
    WORKING_DIRECTORY ${TEST_DIR}
)