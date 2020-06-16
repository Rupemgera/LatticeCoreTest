cmake_minimum_required(VERSION 3.10)

set(project_name CTest)
project(${project_name})

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(googletest)

# if(${CMAKE_VERSION} VERSION_LESS 3.12)
#     file(GLOB TEST_FILES ${PROJECT_SOURCE_DIR}/*.cxx)
# else()
#     file(GLOB TEST_FILES ${PROJECT_SOURCE_DIR}/*.cxx CONFIG_DEPENDS)
# endif()
set(TEST_FILES 
# "mesh_test.cxx"
# "frame_test.cxx"
# "field_test.cxx"
    "FEM_test.cxx"
    "PSO_test.cxx"
)
option(cuda_test_on "wether test CUDA modules" OFF)
if(${cuda_test_on})
    set(${TEST_FILES} ${TEST_FILES} "GPU_test.cxx")
endif()

find_package(dlib REQUIRED)
message("CTets dlib : " ${dlib_FOUND})

add_executable(${project_name} ${TEST_FILES})

enable_testing()

# include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

target_include_directories(${project_name} PUBLIC ${target_path}/LatticeCore)

target_link_libraries(${project_name} PUBLIC  gtest gtest_main dlib::dlib LatticeCore)

add_test(NAME unitTests COMMAND ${PROJECT_NAME})
