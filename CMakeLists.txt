cmake_minimum_required(VERSION 3.10)

project(LatticeCoreTest)

## Eigen
if (NOT Eigen3_FOUND)
	if (UNIX)
		find_package (Eigen3 REQUIRED)
		include_directories(${EIGEN3_INCLUDE_DIR})
		message(${EIGEN3_INCLUDE_DIR})
	endif()
	if(WIN32)
		include_directories("$ENV{EIGEN_DIR}")
		message("Eigen3 directory : " "$ENV{EIGEN_DIR}")
		set(Eigen3_FOUND true)
	endif()
endif()

add_subdirectory(target/LatticeCore)
include_directories("target/LatticeCore")

set(testdata_path "${PROJECT_SOURCE_DIR}/testdata")
add_subdirectory(CTest)