cmake_minimum_required(VERSION 3.30)

# Set the project name and version
project(BandStructure VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find the Armadillo package
find_package(Armadillo REQUIRED)

# Add the executable
add_executable(${PROJECT_NAME} main.cpp)

# Link the Armadillo library
target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIBRARIES})

# Include Armadillo headers
target_include_directories(${PROJECT_NAME} PRIVATE ${ARMADILLO_INCLUDE_DIRS})