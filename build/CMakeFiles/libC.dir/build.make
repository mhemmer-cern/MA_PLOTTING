# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/tavlin/Documents/MA_Plotting

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tavlin/Documents/MA_Plotting/build

# Include any dependencies generated for this target.
include CMakeFiles/libC.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/libC.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/libC.dir/flags.make

# Object files for target libC
libC_OBJECTS =

# External object files for target libC
libC_EXTERNAL_OBJECTS =

liblibC.so: CMakeFiles/libC.dir/build.make
liblibC.so: liblibB.so
liblibC.so: libtestLib.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libCore.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libImt.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libRIO.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libNet.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libHist.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libGraf.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libGraf3d.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libGpad.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libROOTDataFrame.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libTree.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libTreePlayer.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libRint.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libPostscript.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libMatrix.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libPhysics.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libMathCore.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libThread.so
liblibC.so: /home/tavlin/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-02-alice7-7/lib/libMultiProc.so
liblibC.so: CMakeFiles/libC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/tavlin/Documents/MA_Plotting/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX shared library liblibC.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/libC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/libC.dir/build: liblibC.so

.PHONY : CMakeFiles/libC.dir/build

CMakeFiles/libC.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/libC.dir/cmake_clean.cmake
.PHONY : CMakeFiles/libC.dir/clean

CMakeFiles/libC.dir/depend:
	cd /home/tavlin/Documents/MA_Plotting/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tavlin/Documents/MA_Plotting /home/tavlin/Documents/MA_Plotting /home/tavlin/Documents/MA_Plotting/build /home/tavlin/Documents/MA_Plotting/build /home/tavlin/Documents/MA_Plotting/build/CMakeFiles/libC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/libC.dir/depend

