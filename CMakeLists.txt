cmake_minimum_required(VERSION 3.16)

project(Ax LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Create Ax executable
add_executable(Ax
    src/Ax.cpp
)

target_include_directories(Ax
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

# Create Atb executable 
add_executable(Atb
    src/Atb.cpp
)

target_include_directories(Atb
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)