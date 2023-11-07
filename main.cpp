#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{

    // Change the definition here to change resolution
    Scene scene(784, 784);

    Material* red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f);
    Material* green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f);
    Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Material* light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f+0.058f, 0.747f+0.258f, 0.747f) + 15.6f * Vector3f(0.740f+0.287f,0.740f+0.160f,0.740f) + 18.4f *Vector3f(0.737f+0.642f,0.737f+0.159f,0.737f)));
    light->Kd = Vector3f(0.65f);

    Material* whiteM = new Material(MICROFACET, Vector3f(0.0f));
    whiteM->Ks = Vector3f(0.45, 0.45, 0.45);
    whiteM->Kd = Vector3f(0.3, 0.3, 0.25);

    Material* Microfacet = new Material(MICROFACET, Vector3f(0.0f));
    Microfacet->Kd = Vector3f(0.5, 0.5, 0.5);
    Microfacet->Ks = Vector3f(0.5, 0.5, 0.5);
    

    Material* mirror = new Material(MIRROR, Vector3f(0.0f)); 
    mirror->Ks = Vector3f(0.45, 0.45, 0.45); 
    mirror->Kd = Vector3f(0.3, 0.3, 0.25); 
    mirror->ior = 12.85;

    Sphere sphere(Vector3f(420, 90, 130), 90, mirror);

    MeshTriangle floor("C:/Users/32382/Downloads/games08/models/cornellbox/floor.obj", false,white);
    MeshTriangle shortbox("C:/Users/32382/Downloads/games08/models/cornellbox/shortbox.obj", false, white);
    MeshTriangle tallbox("C:/Users/32382/Downloads/games08/models/cornellbox/tallbox.obj", false,mirror);
    MeshTriangle left("C:/Users/32382/Downloads/games08/models/cornellbox/left.obj", false,red);
    MeshTriangle right("C:/Users/32382/Downloads/games08/models/cornellbox/right.obj", false, green);
    MeshTriangle light_("C:/Users/32382/Downloads/games08/models/cornellbox/light.obj", false,light);
    MeshTriangle bunny("C:/Users/32382/Downloads/games08/models/spot/spot_triangulated_good.obj",true, mirror, Vector3f(200,85, 150),Vector3f(140, 140,140));

    scene.Add(&floor);
    scene.Add(&sphere);
    //scene.Add(&shortbox);
    scene.Add(&tallbox);
    scene.Add(&bunny);
    scene.Add(&left);
    scene.Add(&right);
    scene.Add(&light_);

    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}