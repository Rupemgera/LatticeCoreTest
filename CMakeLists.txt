cmake_minimum_required(VERSION 3.10)

project(LatticeCoreTest)

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

# option(target_path "the path to LatticeCore" "target/LatticeCore")
# set(target_path "target/LatticeCore")

if (NOT Eigen3_FOUND)
	if(WIN32)
		include_directories("$ENV{EIGEN_DIR}")
		message("Eigen3 directory : " "$ENV{EIGEN_DIR}")
		set(Eigen3_FOUND true)
	else()
		find_package (Eigen3 REQUIRED)
		# include_directories(${EIGEN3_INCLUDE_DIR})
		message("eigen path: " ${EIGEN3_INCLUDE_DIR})
	endif()
endif()

set(target_path ${CMAKE_SOURCE_DIR}/../LatticeShell/Modules)
message("target " ${target_path})
add_subdirectory(${target_path}/LatticeCore LatticeCore.out)
# include_directories("target/LatticeCore")

add_subdirectory(CTest)

set(testdata_path "${PROJECT_SOURCE_DIR}/testdata")
configure_file(
	"${PROJECT_SOURCE_DIR}/Testdata_path.h.in"
	"${PROJECT_SOURCE_DIR}/CTest/Testdata_path.h"
	NEWLINE_STYLE LF
)
