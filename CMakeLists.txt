cmake_minimum_required(VERSION 3.4)
project(smoothq)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fopenmp -O3 -g")

add_executable(smoothq smoothq.cpp)
