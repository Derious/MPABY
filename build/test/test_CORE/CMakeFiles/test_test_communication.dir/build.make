# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_SOURCE_DIR = /home/zck/文档/MPABY/MPABY

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zck/文档/MPABY/MPABY/build

# Include any dependencies generated for this target.
include test/test_CORE/CMakeFiles/test_test_communication.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/test_CORE/CMakeFiles/test_test_communication.dir/compiler_depend.make

# Include the progress variables for this target.
include test/test_CORE/CMakeFiles/test_test_communication.dir/progress.make

# Include the compile flags for this target's objects.
include test/test_CORE/CMakeFiles/test_test_communication.dir/flags.make

test/test_CORE/CMakeFiles/test_test_communication.dir/test_communication.cpp.o: test/test_CORE/CMakeFiles/test_test_communication.dir/flags.make
test/test_CORE/CMakeFiles/test_test_communication.dir/test_communication.cpp.o: /home/zck/文档/MPABY/MPABY/test/test_CORE/test_communication.cpp
test/test_CORE/CMakeFiles/test_test_communication.dir/test_communication.cpp.o: test/test_CORE/CMakeFiles/test_test_communication.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zck/文档/MPABY/MPABY/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/test_CORE/CMakeFiles/test_test_communication.dir/test_communication.cpp.o"
	cd /home/zck/文档/MPABY/MPABY/build/test/test_CORE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/test_CORE/CMakeFiles/test_test_communication.dir/test_communication.cpp.o -MF CMakeFiles/test_test_communication.dir/test_communication.cpp.o.d -o CMakeFiles/test_test_communication.dir/test_communication.cpp.o -c /home/zck/文档/MPABY/MPABY/test/test_CORE/test_communication.cpp

test/test_CORE/CMakeFiles/test_test_communication.dir/test_communication.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_test_communication.dir/test_communication.cpp.i"
	cd /home/zck/文档/MPABY/MPABY/build/test/test_CORE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zck/文档/MPABY/MPABY/test/test_CORE/test_communication.cpp > CMakeFiles/test_test_communication.dir/test_communication.cpp.i

test/test_CORE/CMakeFiles/test_test_communication.dir/test_communication.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_test_communication.dir/test_communication.cpp.s"
	cd /home/zck/文档/MPABY/MPABY/build/test/test_CORE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zck/文档/MPABY/MPABY/test/test_CORE/test_communication.cpp -o CMakeFiles/test_test_communication.dir/test_communication.cpp.s

# Object files for target test_test_communication
test_test_communication_OBJECTS = \
"CMakeFiles/test_test_communication.dir/test_communication.cpp.o"

# External object files for target test_test_communication
test_test_communication_EXTERNAL_OBJECTS =

bin/test_test_communication: test/test_CORE/CMakeFiles/test_test_communication.dir/test_communication.cpp.o
bin/test_test_communication: test/test_CORE/CMakeFiles/test_test_communication.dir/build.make
bin/test_test_communication: /usr/local/lib/libemp-tool.so
bin/test_test_communication: /usr/lib/x86_64-linux-gnu/libssl.so
bin/test_test_communication: /usr/lib/x86_64-linux-gnu/libcrypto.so
bin/test_test_communication: test/test_CORE/CMakeFiles/test_test_communication.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zck/文档/MPABY/MPABY/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/test_test_communication"
	cd /home/zck/文档/MPABY/MPABY/build/test/test_CORE && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_test_communication.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/test_CORE/CMakeFiles/test_test_communication.dir/build: bin/test_test_communication
.PHONY : test/test_CORE/CMakeFiles/test_test_communication.dir/build

test/test_CORE/CMakeFiles/test_test_communication.dir/clean:
	cd /home/zck/文档/MPABY/MPABY/build/test/test_CORE && $(CMAKE_COMMAND) -P CMakeFiles/test_test_communication.dir/cmake_clean.cmake
.PHONY : test/test_CORE/CMakeFiles/test_test_communication.dir/clean

test/test_CORE/CMakeFiles/test_test_communication.dir/depend:
	cd /home/zck/文档/MPABY/MPABY/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zck/文档/MPABY/MPABY /home/zck/文档/MPABY/MPABY/test/test_CORE /home/zck/文档/MPABY/MPABY/build /home/zck/文档/MPABY/MPABY/build/test/test_CORE /home/zck/文档/MPABY/MPABY/build/test/test_CORE/CMakeFiles/test_test_communication.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/test_CORE/CMakeFiles/test_test_communication.dir/depend

