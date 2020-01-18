cmake_minimum_required(VERSION 3.8)

set(project_name CTest)
project(${project_name})

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_definitions(-DTEST_DATA_PATH="${testdata_path}")

add_subdirectory(googletest)

file(GLOB TEST_FILES ${PROJECT_SOURCE_DIR}/*.cxx)

find_package(dlib REQUIRED)
message("CTets dlib : " ${dlib_FOUND})

add_executable(${project_name} ${TEST_FILES})

enable_testing()

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

target_link_libraries(${project_name} LatticeCore gtest gtest_main dlib::dlib)

add_test(NAME unitTests COMMAND ${PROJECT_NAME})
