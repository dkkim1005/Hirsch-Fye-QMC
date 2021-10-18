# Copyright (c) 2020 Dongkyu Kim (dkkim1005@gmail.com)
CMAKE_MINIMUM_REQUIRED (VERSION 3.10 FATAL_ERROR)

SET (LANGUAGES C CXX)
PROJECT (INSTALL_LIBRARIES LANGUAGES ${LANGUAGES})

# add trng4 library
INCLUDE (${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindTRNG4.cmake)

# add pybind11 library
IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
    SET (PYBIND11_GIT_REPOSITORY https://github.com/pybind/pybind11.git)
    MESSAGE (NOTICE " pybind11 library is cloned at the current source directory.")
    EXECUTE_PROCESS (COMMAND git clone ${PYBIND11_GIT_REPOSITORY} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pybind11/CMakeLists.txt)
        MESSAGE (FATAL_ERROR " We can't find CMakeLists.txt! Please check a git repository again.")
    ENDIF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pybind11/CMakeLists.txt)
ENDIF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
