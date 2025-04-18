# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.20...3.27)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX
  )

# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
set(PYBIND11_PYTHON_VERSION ">=3.10")
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Add a library using FindPython's tooling (pybind11 also provides a helper like
# this)
python_add_library(
  _core
  MODULE
    src/pybind_module.cpp
    src/dubins_pybind.cpp
    src/astar_pybind.cpp
    src/rational_trig_pybind.cpp
    src/dubins_astar_pybind.cpp
  WITH_SOABI
)
target_link_libraries(_core PRIVATE pybind11::headers)

# Add C++ libraries that interface with _core
add_library(dubins_cpp STATIC src/dubins.cpp)
add_library(astar_cpp STATIC src/astar.cpp)
add_library(rational_trig_cpp STATIC src/rational_trig.cpp)
add_library(dubins_astar_cpp STATIC src/dubins_astar.cpp)

# Library linkeages
target_link_libraries(
  dubins_astar_cpp
  PRIVATE
    dubins_cpp
    astar_cpp
    rational_trig_cpp
)
target_link_libraries(
  _core
  PRIVATE
    dubins_cpp
    astar_cpp
    rational_trig_cpp
    dubins_astar_cpp
)

# This is passing in the version as a define just as an example
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

# The install directory is the output (wheel) directory
install(
  TARGETS
    dubins_cpp
    astar_cpp
    rational_trig_cpp
    dubins_astar_cpp
    _core
  DESTINATION
    mpastar
)
