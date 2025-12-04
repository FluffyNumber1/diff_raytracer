# Raytracer

A sophisticated C++ implementation with OpenGL for visualization. This raytracer renders spheres and tetrahedra in 3D space and displays them within the window.
It also supports lighting, reflections, and animation.

## Features
1. Multiple Geometric Primitives: Spheres and Tetrahedrons with triangle-based ray intersection.
2. Advanced Lighting Model:
- Lambertian diffuse lighting
- Blinn-Phong specular highlights
- Ambient lighting
- Shadow casting with occlusion detection
3. Reflective Materials: Configurable glazed surfaces with adjustable reflectivity
## Dependencies

This project requires the following libraries:

- **OpenGL** (system-wide library)
- **GLFW3**
- **GLEW**

## Usage
This project is set up as a **CLion project** with CMake. To run it:
1. Extract the .zip to a desired folder and open it in CLion or IDE with CMake.
2. **Install required libraries**:
    - [GLFW](https://www.glfw.org/download.html)
    - [GLEW](http://glew.sourceforge.net/)
    - OpenGL (already included on most systems)
2.  Update library paths:  
    Inside `CMakeLists.txt`, change the following variables to match the location of **GLEW** and **GLFW** on your system:
   ```cmake
   set(GLEW_DIR "path/to/glew-2.1.0")
   set(GLFW_DIR "path/to/glfw-3.4.bin.WIN64")
    #Also update relevant paths
```
3. Once correct paths for the user's GLFW and GLEW are updated within the `CMakeLists.txt`, the user can simply run:
    ```powershell
    ./raytracer
    ```
    This calls the main function and display the objects within the scene.

## Controls
`P Key`: Once the raytracer is run, toggle between perspective and orthographic projection modes.

## Scene Objects: 
You can modify the parameters within the following cody by modifying the following
code within the `main.cpp`
```powershell
Sphere s1 = Sphere(rgb(1, 0, 0), Vector(0, 5, 2), 1.0, true, 0.1);
```
Edit the rgb values (between 0-1), the position vector, radius, whether the object is glazed or not, and the reflectivity.
```powershell
Tetrahedron t1 = Tetrahedron(
    Vector(-1, 3, -5),  // v0
    Vector(1, 3, -5),   // v1  
    Vector(0, 5, -5),   // v2
    Vector(0, 6, -6),   // v3
    rgb(1, 1, 0),       // Yellow color
    true,               // Glazed (reflective)
    0.4f                // Reflectivity coefficient
);
```
Add a tetrahedron with 4 triangular faces with the above qualities.
```powershell
const Camera camera = Camera(Vector(0, 1, 0), Vector(0,0, 0), Vector(0, 0, 1));
```
You can edit the camera lookAt direction and position.

## Rendering Quality

```powershell
int imageSize = 5000;
```
Higher values = better quality => slower rendering
```powershell
int windowSize = 512; ;
```
Display window resolution.
```powershell
float scale = 3.0f;
```
Adjusts the camera's FOV for greater field and visual range.


## Code Structure:

`rgb` - Color representation class

`Vector` - 3D vector math operations

`Sphere` - Sphere primitive with color, position, and radius

`Triangle` - Triangle primitive for complex geometry

`Tetrahedron` - 4-faced polyhedron composed of triangles

`Camera` - Camera system with position and orientation

## Methods

`getIntersections()` - Ray sphere/triangle intersection calculation

`generateImage()` - OpenGL display and window management

`isInShadow()` - Shadow ray casting for occlusion detection

`rayTriangleIntersectionMullerOptimized()` - Möller-Trumbore algorithm implementation

`reflect` - Mirror reflection calculation

`rayTriangleIntersection()` - standard ray-triangle approach