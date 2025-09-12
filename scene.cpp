#include "scene.h"

void setup_scene(Scene_t& rendererScene) {//ran once to mesh meshes and instance instances
    std::vector<Tri4d_t> mesh1 = {
        // mesh origin (center at 0,0,0)
        {{0,0,0,1}, {0,0,0,1}, {0,0,0,1}},

        // bound centre (also at origin now)
        {{0,0,0,1}, {0,0,0,1}, {0,0,0,1}},

        // Front face (z = +50)
        {{-50,-50, 50,1}, { 50,-50, 50,1}, { 50, 50, 50,1}},
        {{-50,-50, 50,1}, { 50, 50, 50,1}, {-50, 50, 50,1}},

        // Back face (z = -50)
        {{-50,-50,-50,1}, { 50, 50,-50,1}, { 50,-50,-50,1}},
        {{-50,-50,-50,1}, {-50, 50,-50,1}, { 50, 50,-50,1}},

        // Left face (x = -50)
        {{-50,-50,-50,1}, {-50,-50, 50,1}, {-50, 50, 50,1}},
        {{-50,-50,-50,1}, {-50, 50, 50,1}, {-50, 50,-50,1}},

        // Right face (x = +50)
        {{ 50,-50,-50,1}, { 50, 50, 50,1}, { 50,-50, 50,1}},
        {{ 50,-50,-50,1}, { 50, 50,-50,1}, { 50, 50, 50,1}},

        // Top face (y = +50)
        {{-50, 50,-50,1}, {-50, 50, 50,1}, { 50, 50, 50,1}},
        {{-50, 50,-50,1}, { 50, 50, 50,1}, { 50, 50,-50,1}},

        // Bottom face (y = -50)
        {{-50,-50,-50,1}, { 50,-50, 50,1}, {-50,-50, 50,1}},
        {{-50,-50,-50,1}, { 50,-50,-50,1}, { 50,-50, 50,1}},
    };

    struct Mesh_t cubeCube;
    cubeCube.Triangles = mesh1;
    cubeCube.triCount = mesh1.size();
    cubeCube.boundBox = {
        {{-50,-50,-50}},
        {{-50,-50,50}},
        {{-50,50,-50}},
        {{-50,50,50}},
        {{50,-50,-50}},
        {{50,-50,50}},
        {{50,50,-50}},
        {{50,50,50}}
    };

    struct MeshInstance_t mesh1Instance1;
    mesh1Instance1.parentMesh = cubeCube;
    mesh1Instance1.pos = { 0,0,0 };

    struct MeshInstance_t mesh1Instance2;
    mesh1Instance2.parentMesh = cubeCube;
    mesh1Instance2.pos = { 0,50,-150 };

    struct MeshInstance_t mesh1Instance3;
    mesh1Instance3.parentMesh = cubeCube;
    mesh1Instance3.pos = { 100,0,-300 };

    struct Camera_t camera;
    camera.position = { 0,30,300 };//by default facing down -Z
    camera.rotation = { 0,0,0 };

    rendererScene.meshInstances.push_back(mesh1Instance1);
    rendererScene.meshInstances.push_back(mesh1Instance2);
    //rendererScene.meshInstances.push_back(mesh1Instance3);
    rendererScene.meshInstanceCount = rendererScene.meshInstances.size();
    rendererScene.camera = camera;
}