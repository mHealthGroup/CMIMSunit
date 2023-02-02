cmake_minimum_required(VERSION 3.5)
# Project name
project(CMIMS VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(config.h.in config.h)

file(GLOB_RECURSE sources RELATIVE ${CMAKE_SOURCE_DIR} "src/*.c")
list(REMOVE_ITEM sources "src/cli.c")

add_executable(cmims ${sources} "src/cli.c")
add_library(mims STATIC ${sources})

target_include_directories(cmims PUBLIC
                           ${PROJECT_BINARY_DIR} PUBLIC src PUBLIC src/signal PUBLIC src/stats PUBLIC src/trapz
                           )

target_include_directories(mims PUBLIC
                           ${PROJECT_BINARY_DIR} PUBLIC src PUBLIC src/signal PUBLIC src/stats PUBLIC src/trapz
                           )
