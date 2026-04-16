cmake_minimum_required(VERSION 3.10)
project(optimization)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#GCGEN
set(GCGEN_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/GCGen)
set(GCGEN_SAMPLE_INCLUDE ${GCGEN_ROOT}/sample_src)
set(GCGEN_CORE_INCLUDE ${GCGEN_ROOT}/src)
set(GCGEN_SOURCES
    ${GCGEN_ROOT}/sample_src/Hill/HillProblem.cpp
    ${GCGEN_ROOT}/sample_src/Shekel/ShekelProblem.cpp
    ${GCGEN_ROOT}/src/IOptProblem.cpp
    ${GCGEN_ROOT}/src/IGeneralOptProblem.cpp
    ${GCGEN_ROOT}/src/IConstrainedOptProblem.cpp
    ${GCGEN_ROOT}/src/IGeneralOptProblemFamily.cpp
)


# DISLIN
set(DISLIN_INCLUDE_DIR /usr/local/dislin/include)
set(DISLIN_LIB_DIR /usr/local/dislin/lib)

# поиск библиотек
find_library(DISLIN_CPP_LIB discpp PATHS ${DISLIN_LIB_DIR} REQUIRED)
find_library(DISLIN_C_LIB dislin PATHS ${DISLIN_LIB_DIR} REQUIRED)
find_library(X11_LIB X11 REQUIRED)

# основной файл
set(MAIN_TARGET optimization)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${DISLIN_INCLUDE_DIR}
)

link_directories(${DISLIN_LIB_DIR})

add_executable(${MAIN_TARGET} main.cpp)

set_target_properties(${MAIN_TARGET} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${MAIN_TARGET}
    ${DISLIN_CPP_LIB}
    ${DISLIN_C_LIB}
    ${X11_LIB}
    m
    pthread
)

# тестовый файл
set(TEST_TARGET test_gcgen)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${GCGEN_SAMPLE_INCLUDE}
    ${GCGEN_CORE_INCLUDE}
)

add_executable(${TEST_TARGET} test_gcgen.cpp ${GCGEN_SOURCES})

set_target_properties(${TEST_TARGET} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${TEST_TARGET}
    ${DISLIN_CPP_LIB}
    ${DISLIN_C_LIB}
    ${X11_LIB}
    m
    pthread
)

# characteristics
set(CHAR_TARGET characteristics)

add_executable(${CHAR_TARGET} characteristics.cpp ${GCGEN_SOURCES})

set_target_properties(${CHAR_TARGET} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${CHAR_TARGET}
    ${DISLIN_CPP_LIB}
    ${DISLIN_C_LIB}
    ${X11_LIB}
    m
    pthread
)
