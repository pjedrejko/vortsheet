cmake_minimum_required(VERSION 3.14)
project(VortexSheet CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4 /Od)
else()
    add_compile_options(-Wall -Wextra -O3 -march=native)
endif()

#add_compile_options(-fsanitize=address -g)
#add_link_options(-fsanitize=address)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(vortsheet 
    src/main.cpp
    src/Sheet.cpp
    src/Sheet.hpp
    src/Nodes.hpp
    src/Nodes.cpp
    src/Cell.hpp
    src/Cell.cpp
    src/Tree.hpp
    src/Tree.cpp
    src/elliptic.cpp
    src/elliptic.hpp
    src/Velocity.cpp
    src/Velocity.hpp
    src/Params.hpp
    src/Graph.cpp
    src/Graph.hpp
    src/Integrals.cpp
    src/RungeKutta4.cpp
    src/Surgery.cpp
    src/Sync.cpp
)

target_link_libraries(vortsheet PRIVATE OpenMP::OpenMP_CXX MPI::MPI_CXX)
