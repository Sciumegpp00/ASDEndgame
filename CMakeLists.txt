cmake_minimum_required(VERSION 3.19)
project(ASDEndgame)
set(CMAKE_CXX_STANDARD 14)

FILE(GLOB CppSources **/**/*.cpp)
FILE(GLOB HeadersSources **/**/*.h)

add_executable(ASDEndgame ${HeadersSources} ${CppSources} endgame.cpp)