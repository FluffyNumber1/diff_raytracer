cmake_minimum_required(VERSION 3.30)
project(raytracer)

set(CMAKE_CXX_STANDARD 20)

add_executable(raytracer main.cpp)

# GLEW & GLFW directories
#Edit the following directories to the location on your system
set(GLEW_DIR "C:/Users/nadan/glew-2.1.0")
set(GLFW_DIR "C:/Users/nadan/glfw-3.4.bin.WIN64")

# Find OpenGL
find_package(OpenGL REQUIRED)

# Set up include directories (header files)
set(GLFW_INCLUDE_DIR "${GLFW_DIR}/include")
set(GLEW_INCLUDE_DIR "${GLEW_DIR}/include")

# Set up library files
set(GLFW_LIBRARY "${GLFW_DIR}/lib-mingw-w64/libglfw3.a")
set(GLEW_LIBRARY "${GLEW_DIR}/lib/Release/x64/glew32.lib")

# Set up DLL path (for runtime)
set(GLEW_DLL "${GLEW_DIR}/bin/Release/x64/glew32.dll")

include_directories(raytracer
        ${GLFW_INCLUDE_DIR}
        ${GLEW_INCLUDE_DIR}
)

# Link libraries
target_link_libraries(raytracer
        OpenGL::GL
        ${GLFW_LIBRARY}
        ${GLEW_LIBRARY}
)

# Copy GLEW DLL to output directory
add_custom_command(TARGET raytracer POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${GLEW_DLL}"
        $<TARGET_FILE_DIR:raytracer>
)