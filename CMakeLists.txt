cmake_minimum_required(VERSION 3.19)
project(ASDEndgame)

set(CMAKE_CXX_STANDARD 14)

#FILE(GLOB CppSources **/**/*.cpp)
#FILE(GLOB HeadersSources **/**/*.h)

add_executable(endgame endgame.cpp lib/endgame.h lib/grader.cpp)
#add_executable(endgame ${HeadersSources} ${CppSources} endgame.cpp)

#/usr/bin/g++ -DEVAL -std=c++11 -O2 -pipe -static -s -o endgame grader.cpp endgame.cpp