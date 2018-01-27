cmake_minimum_required(VERSION 3.6)
project(pa3_given)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        main.cpp
        png.cpp
        png.h
        quadtree.cpp
        quadtree.h
        quadtree_given.cpp
        quadtree_given.h
        rgbapixel.cpp
        rgbapixel.h)

add_executable(pa3_given ${SOURCE_FILES})