cmake_minimum_required(VERSION 3.4)

# vllib is a shared library so it can be linked to binary extension modules
# that want to use it from C++, rather than from python.

message("SUBDIRECTORY md2")

project(md2 CXX)

set(CMAKE_CXX_STANDARD 17)

add_library( md2
  SHARED
# source file(s)
  verletlist/vl.cpp
#   forces/forces.cpp
  spatial_sorting/spatial_sorting.cpp
  spatial_sorting/hilbert.cpp
  spatial_sorting/hilbert_c.cpp
)

# The shared library is installed at the same location as the binary extension modules 
install(TARGETS md2
    LIBRARY DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/.."
)