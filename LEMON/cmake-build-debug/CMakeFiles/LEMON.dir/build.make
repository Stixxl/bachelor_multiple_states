# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/corbinian/Programs/clion-2018.1.1/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/corbinian/Programs/clion-2018.1.1/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/LEMON.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/LEMON.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LEMON.dir/flags.make

CMakeFiles/LEMON.dir/main.cpp.o: CMakeFiles/LEMON.dir/flags.make
CMakeFiles/LEMON.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/LEMON.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LEMON.dir/main.cpp.o -c /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/main.cpp

CMakeFiles/LEMON.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LEMON.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/main.cpp > CMakeFiles/LEMON.dir/main.cpp.i

CMakeFiles/LEMON.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LEMON.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/main.cpp -o CMakeFiles/LEMON.dir/main.cpp.s

CMakeFiles/LEMON.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/LEMON.dir/main.cpp.o.requires

CMakeFiles/LEMON.dir/main.cpp.o.provides: CMakeFiles/LEMON.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/LEMON.dir/build.make CMakeFiles/LEMON.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/LEMON.dir/main.cpp.o.provides

CMakeFiles/LEMON.dir/main.cpp.o.provides.build: CMakeFiles/LEMON.dir/main.cpp.o


# Object files for target LEMON
LEMON_OBJECTS = \
"CMakeFiles/LEMON.dir/main.cpp.o"

# External object files for target LEMON
LEMON_EXTERNAL_OBJECTS =

LEMON: CMakeFiles/LEMON.dir/main.cpp.o
LEMON: CMakeFiles/LEMON.dir/build.make
LEMON: CMakeFiles/LEMON.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable LEMON"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LEMON.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LEMON.dir/build: LEMON

.PHONY : CMakeFiles/LEMON.dir/build

CMakeFiles/LEMON.dir/requires: CMakeFiles/LEMON.dir/main.cpp.o.requires

.PHONY : CMakeFiles/LEMON.dir/requires

CMakeFiles/LEMON.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LEMON.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LEMON.dir/clean

CMakeFiles/LEMON.dir/depend:
	cd /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/cmake-build-debug /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/cmake-build-debug /home/corbinian/Documents/University/WS17_18/Bachelorthesis/multiple_states/LEMON/cmake-build-debug/CMakeFiles/LEMON.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LEMON.dir/depend
