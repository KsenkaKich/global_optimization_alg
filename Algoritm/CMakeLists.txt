cmake_minimum_required(VERSION 3.10)
project(optimization)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(TARGET_NAME optimization)

set(DISLIN_INCLUDE_DIR /usr/local/dislin/include)
set(DISLIN_LIB_DIR /usr/local/dislin/lib)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${DISLIN_INCLUDE_DIR}
)

link_directories(${DISLIN_LIB_DIR})

find_library(DISLIN_CPP_LIB discpp PATHS ${DISLIN_LIB_DIR} REQUIRED)
find_library(DISLIN_C_LIB dislin PATHS ${DISLIN_LIB_DIR} REQUIRED)
find_library(X11_LIB X11 REQUIRED)

add_executable(${TARGET_NAME} main.cpp)
set_target_properties(${TARGET_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${TARGET_NAME}
    ${DISLIN_CPP_LIB}
    ${DISLIN_C_LIB}
    ${X11_LIB}
    m
    pthread
)