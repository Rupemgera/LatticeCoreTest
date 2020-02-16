
cmake_minimum_required(VERSION 3.10)

set(project_name console)
project(${project_name})

set(src 
"console.cxx"
"eigen.cxx"
)

add_executable(${project_name} ${src})
target_link_libraries(${project_name} PRIVATE LatticeCore)
target_compile_features(LatticeCore PRIVATE cxx_std_14)