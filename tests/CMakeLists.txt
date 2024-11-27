cmake_minimum_required(VERSION 3.27)
add_executable(suite_1_app simulator_test.c)

target_link_libraries(suite_1_app ${LIBRARY_NAME} Unity)

add_test(suite_1_test suite_1_app)