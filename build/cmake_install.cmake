# Install script for directory: /media/tavlin/Samsung_T5/Analyse/MA_Plotting

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELWITHDEBINFO")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/install/plot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/install/plot")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/install/plot"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/install/plot")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/install" TYPE EXECUTABLE FILES "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/build/plot")
  if(EXISTS "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/install/plot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/install/plot")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/install/plot"
         OLD_RPATH "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/build:/home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-08-alice1-6/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/install/plot")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/libtestLib.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/libtestLib.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/libtestLib.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/libtestLib.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib" TYPE SHARED_LIBRARY FILES "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/build/libtestLib.so")
  if(EXISTS "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/libtestLib.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/libtestLib.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/libtestLib.so"
         OLD_RPATH "/home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-08-alice1-6/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/libtestLib.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibB.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibB.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibB.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibB.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib" TYPE SHARED_LIBRARY FILES "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/build/liblibB.so")
  if(EXISTS "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibB.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibB.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibB.so"
         OLD_RPATH "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/build:/home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-08-alice1-6/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibB.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibC.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibC.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibC.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibC.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib" TYPE SHARED_LIBRARY FILES "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/build/liblibC.so")
  if(EXISTS "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibC.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibC.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibC.so"
         OLD_RPATH "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/build:/home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-08-alice1-6/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/media/tavlin/Samsung_T5/Analyse/MA_Plotting/lib/liblibC.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/media/tavlin/Samsung_T5/Analyse/MA_Plotting/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
