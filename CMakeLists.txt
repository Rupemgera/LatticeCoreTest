cmake_minimum_required(VERSION 3.10)

# project(LatticeCoreTest)

# GIT
# find_package(Git QUIET)

# if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
#     execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
#                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                     RESULT_VARIABLE GIT_SUBMOD_RESULT)
#     if(NOT GIT_SUBMOD_RESULT EQUAL "0")
#         message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
#     endif()
# endif()

add_subdirectory(target/LatticeCore)
# include_directories("target/LatticeCore")

add_subdirectory(CTest)
add_subdirectory(console)

set(testdata_path "${PROJECT_SOURCE_DIR}/testdata")
configure_file(
	"${PROJECT_SOURCE_DIR}/Testdata_path.h.in"
	"${PROJECT_SOURCE_DIR}/CTest/Testdata_path.h"
	NEWLINE_STYLE LF
)
configure_file(
	"${PROJECT_SOURCE_DIR}/Testdata_path.h.in"
	"${PROJECT_SOURCE_DIR}/console/Testdata_path.h"
	NEWLINE_STYLE LF
)