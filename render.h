/*
    render.h/render.cpp contains all the functions and constants related to rendering and rasterising a given scene.
    It should only concern itself with turning a scene type into a filled frame buffer
*/
#ifndef RENDER_H
#define RENDER_H

#include "graphics_structs.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif
#include <vector>
#include <algorithm>
#include "SDL3/SDL.h"

extern const int WINDOW_WIDTH;
extern const int WINDOW_HEIGHT;
extern int gFrameBuffer[];
extern double zBuffer[];
extern const float FOV;
extern const double F;
extern const float ASPECT_RATIO;
extern const int NEAR;
extern const int FAR;
extern const float projectionMatrix[4][4];
extern const SDL_PixelFormatDetails* PIXELFORMAT;

/* FUNCTION AND METHOD DECLARATIONS */
void swapVec2d_t(VecP_t* a, VecP_t* b);

int setPixel(int X, int Y, Uint32 Colour);

int setZPixel(int x, int y, double z, Uint32 Colour);

void interpolatei(std::vector<int>& values, int i0, float d0, int i1, float d1);

void interpolatef(std::vector<float>& values, int i0, float d0, int i1, float d1);

void interpolated(std::vector<double>& values, int i0, double d0, int i1, double d1);

Uint32 multiplyColour(Uint32 colour1, float n);

Vec4d_t matrixMult4x4Vec4d_t(const float multMatrix[4][4], Vec4d_t v);

float signedDistance(Plane_t p, Vec4d_t v);

Vec4d_t linePlaneIntersection(Vec4d_t a, Vec4d_t b, Plane_t p);

void projectMeshMatrix(std::vector<Tri4d_t>& Mesh);

std::vector<Tri3d_t> perspectiveDivide(std::vector<Tri4d_t> t);

std::vector<Tri2d_t> viewportTransformation(std::vector<Tri3d_t> t);

void translateMesh(std::vector<Tri4d_t>& Mesh, Vec4d_t Transform);

void rotateMesh(std::vector<Tri4d_t>& Mesh, Vec4d_t Rotation);

void worldSpaceMesh(std::vector<Tri4d_t>& Mesh, MeshInstance_t Instance);

void cameraSpaceMesh(std::vector<Tri4d_t>& Mesh, Camera_t camera);

int clipOrCull(MeshInstance_t mesh, Camera_t camera);

void clipTriangle(std::vector<Tri4d_t>& output, Tri4d_t input, Plane_t plane);

std::vector<Tri4d_t> clipMesh(std::vector<Tri4d_t>& mesh, Camera_t camera);

Vec4d_t crossProduct(Vec4d_t a, Vec4d_t b);

int isFrontFacing(Tri4d_t t);

std::vector<Tri4d_t> backFaceCullMesh(std::vector<Tri4d_t> mesh);

std::vector<Tri4d_t> rotTranClipCullMesh(MeshInstance_t mesh, Camera_t camera);

void drawLine(float x0, float y0, float x1, float y1);

void drawWireTri(Tri2d_t Tri);

void drawFullTri(Tri2d_t Tri, Uint32 Colour);

void drawShadedTri(Tri2d_t Tri, int Colour);

void render_scene(Scene_t scene);

#endif