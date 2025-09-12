/*
    A lot of the structs used here extend beyond the rendering pipeline (such as setting up a scene requiring a scene object)
    if more math related things need to be defined, a graphics_better_name.cpp will be made and populated
*/
#ifndef GRAPHICS_STRUCTS_H
#define GRAPHICS_STRUCTS_H
#include <vector>

/* STRUCTS */
struct Vec4d_t {
    float x = 0;
    float y = 0;
    float z = 0;
    float w = 1;
};

struct Vec3d_t {
    float x = 0;
    float y = 0;
    double z = 0.0;
};

struct VecP_t {//projected vertex in screenspace
    float x = 0;
    float y = 0;
    float h = 1.0;//colour multiplier
    double z = 0.0;//for the z buffer
};

struct Tri4d_t { //we only create these ones
    struct Vec4d_t v0, v1, v2;
};

struct Tri3d_t {
    struct Vec3d_t v0, v1, v2;
};

struct Tri2d_t { //for viewport 3d triangles
    struct VecP_t v0, v1, v2;
};

struct Plane_t {//for clipping planes etc
    Vec4d_t normal;
};

struct Mesh_t {
    std::vector<Tri4d_t> Triangles;
    size_t triCount = 0;
    std::vector<Tri4d_t> boundBox;// 8 vertices but in a triangle due to how stuff is processed here
};

struct MeshInstance_t {//an instance of a predefined mesh to be in scene
    struct Mesh_t parentMesh;
    struct Vec4d_t pos = { 0,0,0 };
    struct Vec4d_t rot = { 0,0,0 };
    struct Vec4d_t scl = { 1,1,1 };
};

struct Camera_t {
    struct Vec4d_t position;
    struct Vec4d_t rotation;
    std::vector<Plane_t> clippingPlanes = {
        {Vec4d_t{1.0, 0.0, 0.0, 1.0}}, //left
        {Vec4d_t{-1.0, 0.0, 0.0, 1.0}}, //right
        {Vec4d_t{0.0, 1.0, 0.0, 1.0}},//bottom
        {Vec4d_t{0.0, -1.0, 0.0, 1.0}},//top 
        {Vec4d_t{0.0, 0.0, -1.0, 1.0}}, //near
        {Vec4d_t{0.0, 0.0, 1.0, 1.0}} //far
    };
};

struct Scene_t {//a scene is a bunch of instances of meshes and a camera
    std::vector<MeshInstance_t> meshInstances;
    size_t meshInstanceCount = 0;
    struct Camera_t camera;
};

#endif