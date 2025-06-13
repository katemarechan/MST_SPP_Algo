cmake_minimum_required(VERSION 3.16)
project(kila)

set(CMAKE_CXX_STANDARD 17)

add_executable(kila
		main.cpp
		MST_SPP.cpp
		Menu.hpp
		ReadFromFile.hpp
)