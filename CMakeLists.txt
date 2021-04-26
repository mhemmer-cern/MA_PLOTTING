cmake_minimum_required(VERSION 3.17)
project(MA_Plotting CXX)
add_executable(plot src/main.cxx)

add_library(testLib SHARED /home/marvin/C_Headers/CommonHeader.h /home/marvin/C_Headers/Plotting_Patrick.h /home/marvin/Documents/git/Header/Plot.h include/plottingheader.hpp)
add_library(libB SHARED include/Systematics.hpp include/PrepareHeader.hpp include/SignalToBackground.hpp include/test.hpp include/YieldScaling.hpp include/SignificanceMC.hpp include/YieldsData.hpp include/PeakComp.hpp include/BeforeScaling.hpp include/FitAfterScaling.hpp include/SameEventToBackgroundRatio.hpp include/Peaks.hpp include/OnlyPeaks.hpp include/ScaleWithUncer.hpp include/Yields.hpp include/Chi2Test.hpp include/QAPlots.hpp)
add_library(libC SHARED include/plotting.hpp)
set_target_properties(testLib PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(libB PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(libC PROPERTIES LINKER_LANGUAGE CXX)

find_package(ROOT)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
  message(STATUS "root  version: ${ROOT_VERSION}")
  target_link_libraries(testLib ${ROOT_LIBRARIES})
else()
  message(ERROR "Could not find root")
endif(ROOT_FOUND)

target_link_libraries(libB testLib)
target_link_libraries(libC libB)
target_link_libraries(plot libC)

install(TARGETS plot
  RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/install
  )
install(TARGETS testLib
  LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
  )

install(TARGETS libB
  LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
  )

install(TARGETS libC
  LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
  )

# by default build optimized code with debug symbols
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()

set(CMAKE_ALLOWED_BUILD_TYPES DEBUG RELEASE RELWITHDEBINFO)
if(NOT CMAKE_BUILD_TYPE IN_LIST CMAKE_ALLOWED_BUILD_TYPES)
  message(FATAL_ERROR "Invalid build type ${CMAKE_BUILD_TYPE}. Use one of: ${CMAKE_ALLOWED_BUILD_TYPES}.")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -DDEBUG -D__DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -ffast-math -DNODEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -ggdb")
message(STATUS "Using CXX flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

#C++ standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "Using C++${CMAKE_CXX_STANDARD}")
