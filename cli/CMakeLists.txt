add_executable(${EXECUTABLE_NAME} "main.c" "cmdline.c")
target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${LIBRARY_NAME})
