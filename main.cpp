#include <GL/glew.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <GLFW/glfw3.h>
using namespace std;

class rgb
{
public:
    float r, g, b;
    rgb(): r(0), g(0), b(0){}
    rgb(float r, float g, float b) : r(r), g(g), b(b) {}
    rgb operator*(float scalar) const
    {
        return rgb(r * scalar, g * scalar, b * scalar);
    }
    rgb operator+(rgb other) const
    {
        return rgb(r + other.r, g + other.g, b + other.b);
    }
};



class Vector
{
public:
    float x, y, z;
    Vector() : x(0), y(0), z(0) {}
    Vector(float x, float y, float z) : x(x), y(y), z(z) {}

    Vector operator+(const Vector& other) const {
        return Vector(x + other.x, y + other.y, z + other.z);
    }
    Vector& operator+=(const Vector& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    Vector operator-(const Vector& other) const {
        return Vector(x - other.x, y - other.y, z - other.z);
    }
    Vector operator*(float scalar) const
    {
        return Vector(x * scalar, y * scalar, z * scalar);
    }

    Vector normalize() const
    {
        float len = sqrt(x*x + y*y + z*z);
        if (len>0.0f)
        {
            return Vector(x/len, y/len, z/len);
        }
        return Vector(0, 0, 0);
    }
};

class Triangle
{
public:
    Vector v0, v1, v2;
    Triangle(const Vector& a, const Vector& b, const Vector& c, const rgb& color)
        : v0(a), v1(b), v2(c) {}
};

class Tetrahedron
{
public:
    rgb color;
    Vector v0, v1, v2, v3;
    vector<Triangle> faces;
    bool glazed;
    float reflectivity;
    Tetrahedron(const Vector& a, const Vector& b, const Vector& c, const Vector& d, const rgb& color, bool glazed, float reflectivity)
       : v0(a), v1(b), v2(c), v3(d), color(color), glazed(glazed), reflectivity(reflectivity)
    {
        faces.push_back(Triangle(v0, v1, v2, color));
        faces.push_back(Triangle(v0, v1, v3, color));
        faces.push_back(Triangle(v0, v2, v3, color));
        faces.push_back(Triangle(v1, v2, v3, color));
    }
};

class Sphere
{
public:
    rgb color;
    Vector center;
    float radius;
    bool glazed;
    float reflectivity;
    Sphere(const rgb& color, const Vector& center, float radius, bool glazed, float reflectivity)
        : color(color), center(center), radius(radius), glazed(glazed), reflectivity(reflectivity) {}
};

Vector cross(const Vector vec1, const Vector vec2)
{
    return Vector(
            vec1.y * vec2.z - vec1.z * vec2.y,
            vec1.z * vec2.x - vec1.x * vec2.z,
            vec1.x * vec2.y - vec1.y * vec2.x);
}

float dot(const Vector vec1, const Vector vec2)
{
    return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
}

class Camera
{
public:
    Vector lookAt;
    Vector viewPoint;
    Vector u, v, w;
    Camera() : lookAt(), viewPoint(), u(), v(), w() {}
    Camera(const Vector& lookAt, const Vector& viewPoint, const Vector& upVector)
        : lookAt(lookAt), viewPoint(viewPoint)
    {

        this->w = (viewPoint-lookAt).normalize();
        this->u = cross(upVector, w).normalize();
        this->v = cross(w,u).normalize();
    }
};

Vector getRayOrigin(const Camera& camera, float u, float v)
{
    Vector rayOrigin = camera.viewPoint + camera.u * u + camera.v * v;
    return rayOrigin;
}
Vector getRayDirection(const Camera& camera, float u, float v)
{
    return (camera.u * u + camera.v * v - camera.w).normalize();
}



void convertToPixels(const vector<vector<rgb>>& image, int width, int height, vector<unsigned char>& pixels) {
    pixels.resize(width * height * 3);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int idx = ((height - 1 - i) * width + j) * 3;
            pixels[idx + 0] = (unsigned char)(std::min(1.0f, image[i][j].r) * 255);
            pixels[idx + 1] = (unsigned char)(std::min(1.0f, image[i][j].g) * 255);
            pixels[idx + 2] = (unsigned char)(std::min(1.0f, image[i][j].b) * 255);
        }
    }
}

const int WIDTH = 64;
const int HEIGHT = 64;

int generateImage(const vector<vector<rgb>>& image)
{
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Raytracer", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return -1;
    }

    // Convert raytraced image to pixels
    std::vector<unsigned char> pixels;
    convertToPixels(image, image[0].size(), image.size(), pixels);


    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        // Draw raytraced image (scaled up to fill the window)
        glPixelZoom((float)WIDTH / image[0].size(), (float)HEIGHT / image.size());
        glDrawPixels(image[0].size(), image.size(), GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
//orig is the starting point of the ray (camera or reflection origin)
//dir is the ray direction (normalized)
//triangle is the triangle with three vertices
//t is the distance along the ray where the intersection occurs
//hitNormal is the surface normal at the intersection point
bool rayTriangleIntersectionMullerOptimized(const Vector& orig, const Vector& dir, const Triangle &tri, float& t, Vector &hitNormal)
{
    const float kEpsilon = 1e-6f; //small threshold for the floating point calculations
    //First build the edges
    Vector v0v1 = tri.v1 - tri.v0;
    Vector v0v2 = tri.v2 - tri.v0;
    //Edge vectors from vertex v0 to v1 and v0 to v2. This will help define the plane of the triangle
    Vector pvec = cross(dir, v0v2);
    //pvec is a helper vector to determine how much v0v1 aligns with the ray
    float det = dot(v0v1, pvec);
    if (fabs(det) < kEpsilon) return false;
    //if det is zero, then there is no intersection
    //if det is negative then that means the triangle is facing the other way relative to the viewing ray
    float invDet = 1.0f / det;
    //inverse of the determinant
    Vector tvec = orig - tri.v0;
    float u = dot(tvec, pvec) * invDet;
    if (u < 0 || u > 1) return false;
    //Logic behind this section is to determine whether the ray actually lands within the triangle.
    //Calculates the barycentric u
    //For a barycentric coord we know that inside the triangle it will sum to 1 and if not inside it will
    Vector qvec = cross(tvec, v0v1);
    float v = dot(dir, qvec) * invDet;
    if (v < 0 || u + v > 1) return false;
    //The step prior checked whether the ray was inside the strip between v0 and v1, now this one checks whether
    //the ray is between v0 and v2 and not past v1 and v2. These two combined will define the inner space of the triangle.
    t = dot(v0v2, qvec) * invDet;
    if (t < kEpsilon) return false;
    hitNormal = cross(v0v1, v0v2).normalize();
    return true;
}
//Simpler method for finding the ray-intersection which uses the triangle-walking method.
bool rayTriangleIntersection(const Vector& orig, const Vector& dir,
                          const Vector& V0, const Vector& V1, const Vector& V2,
                          float &t, Vector &hitNormal)
{
    Vector normal = cross(V1-V0, V2-V0);
    float denom = dot(dir, normal);
    if (abs(denom) <= 1e-6) return false;
    //if the ray direction and the normal are perpendicular, this means that the ray will never intersect the triangle
    t = dot(V0 - orig, normal) / denom;
    if (t < 1e-6f) return false;
    Vector P = orig + dir*t;
    //Edge vectors
    Vector edge0 = V1 - V0;
    Vector edge1 = V2 - V1;
    Vector edge2 = V0 - V2;
    //Vectors from the vertice to the point
    Vector C0 = P - V0;
    Vector C1 = P - V1;
    Vector C2 = P - V2;
    if (dot(normal, cross(edge0,C0))<0) return false;
    if (dot(normal, cross(edge1,C1))<0) return false;
    if (dot(normal, cross(edge2,C2))<0) return false;
    hitNormal = normal.normalize();
    return true;
}



bool isInShadow(const vector<Sphere>& spheres, const vector<Tetrahedron>& tetrahedrons,
                const Vector& hitPoint, const Vector& lightPos)
{
    Vector shadowRayDir = (lightPos - hitPoint).normalize();
    float distanceToLight = sqrt(dot(lightPos - hitPoint, lightPos - hitPoint));

    // Check spheres
    for (const Sphere& sphere : spheres)
    {
        Vector p = hitPoint - sphere.center;
        float a = dot(shadowRayDir, shadowRayDir);
        float b = 2.0f * dot(p, shadowRayDir);
        float c = dot(p, p) - sphere.radius * sphere.radius;
        float discriminant = b*b - 4*a*c;
        if (discriminant >= 0)
        {
            float sqrtDisc = sqrt(discriminant);
            float t1 = (-b - sqrtDisc) / (2*a);
            float t2 = (-b + sqrtDisc) / (2*a);
            float t = (t1 > 0.001f) ? t1 : ((t2 > 0.001f) ? t2 : -1);
            if (t > 0.001f && t < distanceToLight)
                return true;
        }
    }

    // Check tetrahedrons (triangles)
    for (const Tetrahedron& tetra : tetrahedrons)
    {
        for (const Triangle& tri : tetra.faces)
        {
            float t;
            Vector normal;
            if (rayTriangleIntersection(hitPoint, shadowRayDir, tri.v0, tri.v1, tri.v2, t, normal))
            {
                if (t > 0.001f && t < distanceToLight)
                    return true;
            }
        }
    }

    return false;
}
Vector reflect(Vector V, Vector normal)
{
    return (V - normal * 2 * dot(V, normal)).normalize();
}

pair<bool,rgb> getIntersections(const vector<Sphere>& spheres, const vector<Tetrahedron>& tetrahedrons, Vector rayOrigin, Vector rayDirection, Vector lightPos, int depth = 0)
{
    if (depth > 5)
    {
        return make_pair(false, rgb(0,0,0));
    }

    bool intersects = false;
    float closestT = INFINITY;
    rgb hitColor = {0,0,0};
    const Sphere* hitSphere = nullptr;
    const Tetrahedron* hitTetrahedron = nullptr;
    Vector hitPoint, hitNormal;

    // Check sphere intersections
    for (const Sphere& sphere : spheres)
    {
        Vector p = rayOrigin - sphere.center;

        float a = dot(rayDirection, rayDirection);
        float b = 2.0f * dot(p, rayDirection);
        float c = dot(p, p) - sphere.radius * sphere.radius;

        float discriminant = b * b - 4 * a * c;

        if (discriminant >= 0)
        {
            float sqrtDisc = sqrt(discriminant);
            float t1 = (-b - sqrtDisc) / (2 * a);
            float t2 = (-b + sqrtDisc) / (2 * a);

            float t = (t1 > 0.001f) ? t1 : ((t2 > 0.001f) ? t2 : -1);

            if (t > 0.001f && t < closestT)
            {
                closestT = t;
                hitSphere = &sphere;
                hitTetrahedron = nullptr; // Reset tetrahedron hit
                hitPoint = rayOrigin + rayDirection * t;
                hitNormal = (hitPoint - sphere.center).normalize();
                intersects = true;
            }
        }
    }

    // Check tetrahedron intersections
    for (const Tetrahedron& tetra : tetrahedrons)
    {
        // Check intersection with each face of the tetrahedron
        for (const Triangle& face : tetra.faces)
        {
            float t;
            Vector normal;

            if (rayTriangleIntersection(rayOrigin, rayDirection, face.v0,face.v1,face.v2, t, normal))
            {
                if (t > 0.001f && t < closestT)
                {
                    closestT = t;
                    hitSphere = nullptr; // Reset sphere hit
                    hitTetrahedron = &tetra;
                    hitPoint = rayOrigin + rayDirection * t;
                    hitNormal = normal;
                    intersects = true;
                }
            }
        }
    }

    if (intersects)
    {
        Vector lightDirection = (lightPos - hitPoint).normalize();

        // Lambertian effect
        float Ld = max(0.0f, dot(hitNormal, lightDirection));
        // Ambient Constant
        float La = 0.1f;

        // Specular Lighting
        Vector vE = (rayOrigin - hitPoint).normalize(); // Eye vector
        Vector vL = (lightPos - hitPoint).normalize();   // Light vector
        Vector vH = (vL + vE).normalize();
        float Ls = pow(max(0.0f, dot(hitNormal, vH)), 60.0f);

        // Determine the color and material properties based on what was hit
        rgb objectColor;
        bool isGlazed;
        float reflectivity;

        if (hitSphere != nullptr)
        {
            objectColor = hitSphere->color;
            isGlazed = hitSphere->glazed;
            reflectivity = hitSphere->reflectivity;
        }
        else if (hitTetrahedron != nullptr)
        {
            objectColor = hitTetrahedron->color;
            isGlazed = hitTetrahedron->glazed;
            reflectivity = hitTetrahedron->reflectivity;
        }

        if (!isInShadow(spheres,tetrahedrons, hitPoint, lightPos))
        {
            hitColor = objectColor * (La + Ld + Ls);
        }
        else // When in shadow, only show ambient lighting.
        {
            hitColor = objectColor * La;
        }

        if (isGlazed && depth < 5)
        {
            Vector reflectedRay = reflect(rayDirection, hitNormal);
            // Move slightly away from surface to avoid self-intersection
            Vector reflectionOrigin = hitPoint + hitNormal * 0.001f;

            pair<bool, rgb> reflectionResult = getIntersections(spheres, tetrahedrons, reflectionOrigin, reflectedRay, lightPos, depth + 1);

            if (reflectionResult.first)
            {
                hitColor = hitColor * (1.0f - reflectivity) + reflectionResult.second * reflectivity;
            }
        }
    }

    return make_pair(intersects, hitColor);
}
bool perspective = true;
bool front = true;
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_P && action == GLFW_PRESS)
    {
        cout<< "Perspective changed" << endl;
        perspective = !perspective;
    }
}


int main()
{

    Sphere s1 = Sphere(rgb(1, 0, 0), Vector(0, 5, 0), 1.0, true, 0.1);  // Red sphere
    // Sphere s2 = Sphere(rgb(0, 1, 0), Vector(2, 9, 0), 2.0, true, 0.3);   // Green sphere (static)
    // Sphere s3 = Sphere(rgb(0, 0, 1), Vector(-2, 25, 0), 1.5, true, .2);  // Blue sphere (static)

    vector<Tetrahedron> tetrahedrons;
    Tetrahedron t1 = Tetrahedron(
    Vector(-1, 3, -5),  // v0
    Vector(1, 3, -5),   // v1
    Vector(0, 5, -5),   // v2
    Vector(0, 6, -6),   // v3
    rgb(1, 1, 0),       // Yellow
    true,
    0.4f
);
    // tetrahedrons.push_back(t1);

    Tetrahedron t2 = Tetrahedron(
        Vector(2, 1, -7),   // v0
        Vector(4, 1, -7),   // v1
        Vector(3, 3, -7),   // v2
        Vector(3, 4, -5.5), // v3
        rgb(0, 1, 1),       // Cyan
        true,
        0.3f
    );
    // tetrahedrons.push_back(t2);
    Camera camera = Camera(Vector(0, 1, 0), Vector(0,0, 0), Vector(0, 0, 1)); //Angle towards spheres
    //Camera camera = Camera(Vector(0,3,-5), Vector(0,3,0), Vector(0,1,0));



    int windowSize = 512;
    int imageSize = 128;
    vector<vector<rgb>> image(imageSize, std::vector<rgb>(imageSize));

    if (!glfwInit())
    {
        cerr << "Failed to initialize GLFW" << endl;
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(windowSize, windowSize, "Animated Raytracer", NULL, NULL);
    if (!window) {
        cerr << "Failed to create window" << endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    if (glewInit() != GLEW_OK) {
        cerr << "Failed to initialize GLEW" << endl;
        return -1;
    }

    // float time = 0.0f;
    // float deltaTime = 0.1f;

    Vector lightCenter = Vector(10, 0, 0);
    float lightRadius = 8.0f;


    Vector sphereCenter = Vector(0, 5, -2);
    float sphereRadius = 2.0f;

    while (!glfwWindowShouldClose(window))
    {
        // Update time
        // time += deltaTime;


        Vector lightPos = lightCenter;

        // Vector lightPos = Vector(
        //     lightCenter.x + lightRadius * cos(time),
        //     lightCenter.y + lightRadius * sin(time),
        //     lightCenter.z
        // );
        //
        // Vector redSpherePos = Vector(
        //     sphereCenter.x + sphereRadius * cos(time * 0.7f),
        //     sphereCenter.y,
        //     sphereCenter.z + sphereRadius * sin(time * 0.7f)
        // );
        glfwSetKeyCallback(window, key_callback);

        // s1.center = redSpherePos;

        vector<Sphere> objects = {s1};

        for (int i = 0; i < imageSize; i++)
        {
            for (int j = 0; j < imageSize; j++)
            {
                float scale = 1.0f;
                if (!perspective)
                {
                    scale = 2.5f;
                }
                float u = ((2.0f *(.5+ j) / (imageSize)) - 1.0f)*scale;
                float v = ((2.0f *(.5+ i) / (imageSize)) - 1.0f)*scale;

                //Orthogonal projection
                Vector rayOrigin = getRayOrigin(camera, u, v);
                Vector rayDirection = (camera.lookAt).normalize();

                if (perspective)
                {
                    rayOrigin = camera.viewPoint;
                    rayDirection = getRayDirection(camera,u,v);
                }

                pair<bool, rgb> result = getIntersections(objects,tetrahedrons, rayOrigin, rayDirection, lightPos);
                image[i][j] = result.second;
            }
        }

        // Convert and display the image
        vector<unsigned char> pixels;
        convertToPixels(image, imageSize, imageSize, pixels);

        glClear(GL_COLOR_BUFFER_BIT);
        glPixelZoom((float)windowSize / imageSize, (float)windowSize / imageSize);
        glDrawPixels(imageSize, imageSize, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

