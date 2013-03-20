# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mcifti/VolumeRendering/VolRender

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mcifti/VolumeRendering/VolRender

# Include any dependencies generated for this target.
include CMakeFiles/VolRender.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/VolRender.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/VolRender.dir/flags.make

CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o: CMakeFiles/VolRender.dir/flags.make
CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o: src/Raw3DData.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/mcifti/VolumeRendering/VolRender/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o -c /home/mcifti/VolumeRendering/VolRender/src/Raw3DData.cpp

CMakeFiles/VolRender.dir/src/Raw3DData.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VolRender.dir/src/Raw3DData.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/mcifti/VolumeRendering/VolRender/src/Raw3DData.cpp > CMakeFiles/VolRender.dir/src/Raw3DData.cpp.i

CMakeFiles/VolRender.dir/src/Raw3DData.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VolRender.dir/src/Raw3DData.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/mcifti/VolumeRendering/VolRender/src/Raw3DData.cpp -o CMakeFiles/VolRender.dir/src/Raw3DData.cpp.s

CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o.requires:
.PHONY : CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o.requires

CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o.provides: CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o.requires
	$(MAKE) -f CMakeFiles/VolRender.dir/build.make CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o.provides.build
.PHONY : CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o.provides

CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o.provides.build: CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o

CMakeFiles/VolRender.dir/src/VolRender.cpp.o: CMakeFiles/VolRender.dir/flags.make
CMakeFiles/VolRender.dir/src/VolRender.cpp.o: src/VolRender.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/mcifti/VolumeRendering/VolRender/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/VolRender.dir/src/VolRender.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/VolRender.dir/src/VolRender.cpp.o -c /home/mcifti/VolumeRendering/VolRender/src/VolRender.cpp

CMakeFiles/VolRender.dir/src/VolRender.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VolRender.dir/src/VolRender.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/mcifti/VolumeRendering/VolRender/src/VolRender.cpp > CMakeFiles/VolRender.dir/src/VolRender.cpp.i

CMakeFiles/VolRender.dir/src/VolRender.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VolRender.dir/src/VolRender.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/mcifti/VolumeRendering/VolRender/src/VolRender.cpp -o CMakeFiles/VolRender.dir/src/VolRender.cpp.s

CMakeFiles/VolRender.dir/src/VolRender.cpp.o.requires:
.PHONY : CMakeFiles/VolRender.dir/src/VolRender.cpp.o.requires

CMakeFiles/VolRender.dir/src/VolRender.cpp.o.provides: CMakeFiles/VolRender.dir/src/VolRender.cpp.o.requires
	$(MAKE) -f CMakeFiles/VolRender.dir/build.make CMakeFiles/VolRender.dir/src/VolRender.cpp.o.provides.build
.PHONY : CMakeFiles/VolRender.dir/src/VolRender.cpp.o.provides

CMakeFiles/VolRender.dir/src/VolRender.cpp.o.provides.build: CMakeFiles/VolRender.dir/src/VolRender.cpp.o

# Object files for target VolRender
VolRender_OBJECTS = \
"CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o" \
"CMakeFiles/VolRender.dir/src/VolRender.cpp.o"

# External object files for target VolRender
VolRender_EXTERNAL_OBJECTS =

VolRender: CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o
VolRender: CMakeFiles/VolRender.dir/src/VolRender.cpp.o
VolRender: CMakeFiles/VolRender.dir/build.make
VolRender: /usr/lib/x86_64-linux-gnu/libQtGui.so
VolRender: /usr/lib/x86_64-linux-gnu/libQtCore.so
VolRender: /usr/lib/x86_64-linux-gnu/libX11.so
VolRender: /usr/lib/x86_64-linux-gnu/libglut.so
VolRender: /usr/lib/x86_64-linux-gnu/libXmu.so
VolRender: /usr/lib/x86_64-linux-gnu/libXi.so
VolRender: /usr/lib/x86_64-linux-gnu/libGLEW.so
VolRender: CMakeFiles/VolRender.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable VolRender"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/VolRender.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/VolRender.dir/build: VolRender
.PHONY : CMakeFiles/VolRender.dir/build

CMakeFiles/VolRender.dir/requires: CMakeFiles/VolRender.dir/src/Raw3DData.cpp.o.requires
CMakeFiles/VolRender.dir/requires: CMakeFiles/VolRender.dir/src/VolRender.cpp.o.requires
.PHONY : CMakeFiles/VolRender.dir/requires

CMakeFiles/VolRender.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/VolRender.dir/cmake_clean.cmake
.PHONY : CMakeFiles/VolRender.dir/clean

CMakeFiles/VolRender.dir/depend:
	cd /home/mcifti/VolumeRendering/VolRender && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mcifti/VolumeRendering/VolRender /home/mcifti/VolumeRendering/VolRender /home/mcifti/VolumeRendering/VolRender /home/mcifti/VolumeRendering/VolRender /home/mcifti/VolumeRendering/VolRender/CMakeFiles/VolRender.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/VolRender.dir/depend

