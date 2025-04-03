# CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

# Project definition
project(colloidgel
    VERSION 1.0.0
    DESCRIPTION "Colloidal gel simulation program"
    LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type configuration
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -O3)
endif()

# Find required packages
find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

# Define source files
file(GLOB_RECURSE SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# Define header files
file(GLOB_RECURSE HEADERS 
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)

# Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Include directories
target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Link libraries
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        GSL::gsl
        GSL::gslcblas
        Boost::program_options
)

# Installation rules
include(GNUInstallDirs)

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES ${HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

# # Testing configuration
# option(BUILD_TESTING "Build tests" ON)
# if(BUILD_TESTING)
#     enable_testing()
#     add_subdirectory(tests)
# endif()

# Documentation configuration
option(BUILD_DOCS "Build documentation" OFF)
if(BUILD_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        add_subdirectory(docs)
    endif()
endif()

# Export targets
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)