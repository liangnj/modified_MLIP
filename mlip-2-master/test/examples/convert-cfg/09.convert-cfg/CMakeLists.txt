set(TEST_NAME "09.convert-cfg")
set(TEST_DIR ${TEST_BASE_DIR}/examples/${TEST_NAME})
file(MAKE_DIRECTORY ${TEST_DIR}/out)
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/OUTCAR
  ${CMAKE_CURRENT_SOURCE_DIR}/correct_POSCAR
  )
file(COPY ${TEST_SOURCES} DESTINATION ${TEST_DIR})

add_test(
  NAME ${TEST_NAME}
  COMMAND sh -c "($<TARGET_FILE:mlp> convert-cfg OUTCAR ./out/POSCAR --input-format=vasp-outcar --output-format=vasp-poscar > /dev/null; diff correct_POSCAR ./out/POSCAR) > ./out/stdout.log 2>&1"
  WORKING_DIRECTORY ${TEST_DIR}
)