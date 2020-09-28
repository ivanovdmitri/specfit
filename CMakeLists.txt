#
#  D. Ivanov <dmiivanov@gmail.com>
#  09/2020
#
# For building specfit using cmake:
# mkdir obj
# cd obj
# cmake ..
# make -j3
#
cmake_minimum_required(VERSION 3.9)
project(SPECFIT)

# so that the build commands are always shown
set(CMAKE_VERBOSE_MAKEFILE on)

### where to put the include files
set(INCLUDE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/inc")

### to pu the source files (needed for HTML doc)
set(SOURCE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/src")

### where to put the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_PATH}")
set(CMAKE_INSTALL_LIBDIR "${LIBRARY_OUTPUT_PATH}")

### where to put the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_BINDIR "${EXECUTABLE_OUTPUT_PATH}")
set(CMAKE_BINARY_OUTPUT_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

### add CERN ROOT
find_package(ROOT 5.34.01 REQUIRED)
include(${ROOT_USE_FILE})

# don't want GCC complaints about overloaded virtual methods 
# which happens with older versions of ROOT
string(REGEX REPLACE "\\." "" ROOT_INT_VERSION ${ROOT_VERSION})
If(ROOT_INT_VERSION LESS 60000)
string (REPLACE "-Woverloaded-virtual" ""     CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
Endif(ROOT_INT_VERSION LESS 60000)

# add extra warnings flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")

# set the include directories
include_directories(${INCLUDE_OUTPUT_PATH})


# install C++ header files, LinkDef file is treated sperately
# because it's only needed for the ROOT dictionary generation
file(GLOB SPECFIT_HEADERS inc/*.h)
list(FILTER SPECFIT_HEADERS EXCLUDE REGEX ".*LinkDef.h$")
file(INSTALL ${SPECFIT_HEADERS} DESTINATION ${INCLUDE_OUTPUT_PATH})

# build the C++ shared library with ROOT dictionaries
ROOT_GENERATE_DICTIONARY(specfitDict ${CMAKE_SOURCE_DIR}/inc/specfit.h 
  LINKDEF ${CMAKE_SOURCE_DIR}/inc/specfitLinkDef.h)

# listing all SPECFIT source files here
set(SPECFIT_SOURCES
  src/specfit_canv.cxx
  src/specfit_uti.cxx
  src/TBPLF1.cxx
  src/TCRFlux.cxx
  src/TCRFluxFit.cxx
  src/TSPECFITF1.cxx)

# needed for being able to generate full HTML documentation in the build directory 
file(INSTALL ${SPECFIT_SOURCES} DESTINATION ${SOURCE_OUTPUT_PATH})

# function to prefix files; late versions of cmake have this function
# built it but this shuold work with cmake 3.9 or higher. 
FUNCTION(PREPEND var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

# the same source files are used for compilation 
# as those that will be in the HTML documentation 
PREPEND(SPECFIT_INSTALLED_SOURCES "${CMAKE_BINARY_DIR}/" ${SPECFIT_SOURCES})

# the spctrum fit library
add_library(specfit SHARED 
  ${SPECFIT_INSTALLED_SOURCES}
  specfitDict.cxx)
target_link_libraries(specfit ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lMinuit)

# Python library files
set(PYFILES
  encorr_functions.py
  flux_functions.py
  specfit_cpplib.py
  )
file(INSTALL ${PYFILES} DESTINATION ${CMAKE_BINARY_DIR})

# Python executable files
set(PYEXEFILES specfit.py )
file (INSTALL ${PYEXEFILES} 
  DESTINATION ${CMAKE_BINARY_DIR}
  FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
)

# ROOT macros
set(ROOT_MACRO_FILES load_specfit_cpplib.C generate_specfit_htmldoc.C)
file(INSTALL ${ROOT_MACRO_FILES} DESTINATION ${CMAKE_BINARY_DIR})

# HTML documentation
add_custom_target(htmldoc ALL DEPENDS specfit ${ROOT_MACRO_FILES})
add_custom_command(
    TARGET htmldoc
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${ROOT_EXECUTABLE} -l -b -q generate_specfit_htmldoc.C >& /dev/null
)


# environmental file
file(INSTALL specfit_env.sh DESTINATION ${CMAKE_BINARY_DIR})

# simulations directory
file(INSTALL sim DESTINATION ${CMAKE_BINARY_DIR})

# data directory
file(INSTALL data DESTINATION ${CMAKE_BINARY_DIR})
