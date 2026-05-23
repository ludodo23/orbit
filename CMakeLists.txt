cmake_minimum_required(VERSION 3.16)

file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" PROJECT_VERSION)

project(orbit VERSION ${PROJECT_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -------------------------
# LIB HEADER-ONLY
# -------------------------
add_library(orbit INTERFACE)
target_include_directories(orbit INTERFACE include/)

# -------------------------
# TESTS
# -------------------------
enable_testing()
add_subdirectory(tests)