#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif
#include "SDL3/SDL.h"
#include "SDL3/SDL_main.h"

SDL_Window* gSDLWindow;
SDL_Renderer* gSDLRenderer; 
SDL_Texture* gSDLTexture;
static int gDone;
const int WINDOW_WIDTH = 1920;
const int WINDOW_HEIGHT = 1080;
int* gFrameBuffer;
double* zBuffer;
const float FOV = 90.0f * 3.14159f / 180.0f;//60 degree fov
const double F = 1.0f / tan(FOV / 2.0f);
const float ASPECT_RATIO = (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT;
const int NEAR = 10;
const int FAR = 2000;
const float projectionMatrix[4][4] = {
    {F / ASPECT_RATIO, 0, 0, 0},
    {0, F, 0, 0},
    {0, 0, -((FAR + NEAR) / (NEAR - FAR)), -((2 * FAR * NEAR) / (NEAR - FAR))},
    {0, 0, -1, 0}
};
const SDL_PixelFormatDetails* PIXELFORMAT = SDL_GetPixelFormatDetails(SDL_PIXELFORMAT_RGBA8888);


/*
TODO
19-30fps
53-84fps
use bitshifts for colour mult
better early rejection for h_segment where its all 1.0
hidden surface removal
occulision culling
shading
textures
*/

struct Vec4d_t {
    float x, y, z = 0.0;
    float w = 1;
};

struct Vec3d_t{
    float x, y = 0.0;
    double z = 0.0;
}; 

struct VecP_t {//projected vertex in screenspace
    float x, y = 0.0;
    float h = 1.0;//colour multiplier
    double z = 0.0;//for the z buffer
};

struct Tri4d_t{ //we only create these ones
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
    int triCount;
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
    int meshInstanceCount = 0;
    struct Camera_t camera;
};

/*global variable scene*/
Scene_t rendererScene;

void swapVec2d_t(VecP_t* a, VecP_t* b) {//in go pointer to passed thingy so it doesnt need wierd returning, swap(&a, &b); )
    VecP_t temp = *a;
    *a = *b;
    *b = temp;
}

int setPixel(int X, int Y, Uint32 Colour) {//sets a pixel in the frame buffer to a colour, w/ some limits to avoid index errors
    if (Y == 0) Y++; if (X == 0) X++;
    if (X < 0 || X > WINDOW_WIDTH - 1)
        return -1;
    if (Y < 0 || Y > WINDOW_HEIGHT - 1)
        return -1;
    
    gFrameBuffer[X + ((Y - 1) * WINDOW_WIDTH)] = Colour;
    return 1;
}

int setZPixel(int x, int y, double z, Uint32 Colour) {//z will be inverted as 1/z
    if (y == 0) y++; if (x == 0) x++;
    int pos = x + ((y - 1) * WINDOW_WIDTH);
    if (z > (zBuffer[pos]) + FLT_EPSILON) {
        setPixel(x, y, Colour);
        zBuffer[pos] = z;
        return 1;
    }
    return 0;
}

void interpolatei(std::vector<int> &values, int i0, float d0, int i1, float d1){//pass in a vector interpolates whole integers
    if (i0 == i1) {
        values.clear();
        values.push_back(d0);
    }
    size_t size = abs(i1 - i0) + 1;
    values.resize(size);
    float a = (d1 - d0) / (i1 - i0);
    float d = d0;
    int i, j;
    for (i = i0, j = 0; i <= (int) i1; i++, j++) {
        values[j] = d;
        d += a;
    }
}

void interpolatef(std::vector<float> &values, int i0, float d0, int i1, float d1) {//pass in a vector interpolates floats
    if (i0 == i1) {
        values.clear();
        values.push_back(d0);
    }
    size_t size = abs(i1 - i0) + 1;
    values.resize(size);
    float a = (d1 - d0) / (i1 - i0);
    float d = d0;
    int i, j;
    for (i = i0, j = 0; i <= (int) i1; i++, j++) {
        values[j] = d;
        d += a;
    }
}

void interpolated(std::vector<double>& values, int i0, double d0, int i1, double d1) {//pass in a vector interpolates doubles
    if (i0 == i1) {
        values.clear();
        values.push_back(d0);
    }
    size_t size = abs(i1 - i0) + 1;
    values.resize(size);
    double a = (d1 - d0) / (i1 - i0);
    double d = d0;
    int i, j;
    for (i = i0, j = 0; i <= (int) i1; i++, j++) {
        values[j] = d;
        d += a;
    }
}

Uint32 multiplyColour(Uint32 colour1, float n) {//multiplies two colours based off RGBA format
    if (n == 1) return colour1;
    Uint8 r1, g1, b1, a1;
    SDL_GetRGBA(colour1, PIXELFORMAT, NULL, &r1, &g1, &b1, &a1);
    Uint8 r = (r1 * n);
    Uint8 g = (g1 * n);
    Uint8 b = (b1 * n);
    Uint8 a = (a1 * n);
    return(SDL_MapRGBA(PIXELFORMAT, NULL, r, g, b, a));
}

Vec4d_t matrixMult4x4Vec4d_t(const float multMatrix[4][4], Vec4d_t v) {//matrix multiplication of a 4x4 matrix and a 1x4 coordinate
    int i, j;
    Vec4d_t tv = { 0,0,0,0 };
    for (i = 0; i < 4; i++) {
        float sum = 0;
        for (j = 0; j < 4; j++) {
            sum += *(float*)((unsigned char*)&v + (j * 4)) * multMatrix[i][j];
        }
        *(float*)((unsigned char*)&tv + (i * 4)) = sum;
    }
    
    return tv;
}

float signedDistance(Plane_t p, Vec4d_t v) {//signed distance between a point and a plane in 4d space
    float distance = (v.x * p.normal.x) + (v.y * p.normal.y) + (v.z * p.normal.z) + (v.w * p.normal.w);
    return distance;
}

Vec4d_t linePlaneIntersection(Vec4d_t a, Vec4d_t b, Plane_t p) {//intersection of a plane and two points in 4d space
    float t, sum1, sum2;
    Vec4d_t i;
    sum1 = (a.x * p.normal.x) + (a.y * p.normal.y) + (a.z * p.normal.z) + (a.w * p.normal.w);
    sum2 = (b.x * p.normal.x) + (b.y * p.normal.y) + (b.z * p.normal.z) + (b.w * p.normal.w);
    if (sum2 == 0) {//the line is parallel with the plane
        return a;
    }
    t = sum1 / (sum1 - sum2);
    i.x = a.x + (b.x - a.x) * t;
    i.y = a.y + (b.y - a.y) * t;
    i.z = a.z + (b.z - a.z) * t;
    i.w = a.w + (b.w - a.w) * t;
    return i;
}

void projectMeshMatrix(std::vector<Tri4d_t>& Mesh){//projects the matrix according to the projection matrix...
    
    int i;
    for (i = 0; i < Mesh.size(); i++) {
        Mesh[i].v0 = matrixMult4x4Vec4d_t(projectionMatrix, Mesh[i].v0);
        Mesh[i].v1 = matrixMult4x4Vec4d_t(projectionMatrix, Mesh[i].v1);
        Mesh[i].v2 = matrixMult4x4Vec4d_t(projectionMatrix, Mesh[i].v2);
    }
}

std::vector<Tri3d_t> perspectiveDivide(std::vector<Tri4d_t> t) {//scales the coordinates in accordance to w
    std::vector<Tri3d_t> output;
    output.resize(t.size());
    int i;
    for (i = 0; i < t.size(); i++) {
        output[i].v0.x = t[i].v0.x / t[i].v0.w;
        output[i].v0.y = t[i].v0.y / t[i].v0.w;
        output[i].v0.z = t[i].v0.z / t[i].v0.w;

        output[i].v1.x = t[i].v1.x / t[i].v1.w;
        output[i].v1.y = t[i].v1.y / t[i].v1.w;
        output[i].v1.z = t[i].v1.z / t[i].v1.w;

        output[i].v2.x = t[i].v2.x / t[i].v2.w;
        output[i].v2.y = t[i].v2.y / t[i].v2.w;
        output[i].v2.z = t[i].v2.z / t[i].v2.w;

    }
    return output;
}

std::vector<Tri2d_t> viewportTransformation(std::vector<Tri3d_t> t) {//scales coordinates based on viewport size
    std::vector<Tri2d_t> output;
    output.resize(t.size());
    int i;
    double sum = 0;
    for (i = 0; i < t.size(); i++) {
        output[i].v0.x = (t[i].v0.x + 1) / 2 * WINDOW_WIDTH; 
        output[i].v0.y = (1 - t[i].v0.y) / 2 * WINDOW_HEIGHT;
        output[i].v1.x = (t[i].v1.x + 1) / 2 * WINDOW_WIDTH; 
        output[i].v1.y = (1 - t[i].v1.y) / 2 * WINDOW_HEIGHT;
        output[i].v2.x = (t[i].v2.x + 1) / 2 * WINDOW_WIDTH;
        output[i].v2.y = (1 - t[i].v2.y) / 2 * WINDOW_HEIGHT;

        output[i].v0.z = (1 / -t[i].v0.z);
        output[i].v1.z = (1 / -t[i].v1.z);
        output[i].v2.z = (1 / -t[i].v2.z);
        sum += output[i].v0.z;
        sum += output[i].v1.z;
        sum += output[i].v2.z;
    }
    sum = sum / (t.size() * 3);
    return output;
}

void translateMesh(std::vector<Tri4d_t>& Mesh, Vec4d_t Transform) {//adds two matricies, for translation
    int i;
    for (i = 0; i < Mesh.size(); i++) {
        Mesh[i].v0.x += Transform.x;
        Mesh[i].v0.y += Transform.y;
        Mesh[i].v0.z += Transform.z;

        Mesh[i].v1.x += Transform.x;
        Mesh[i].v1.y += Transform.y;
        Mesh[i].v1.z += Transform.z;

        Mesh[i].v2.x += Transform.x;
        Mesh[i].v2.y += Transform.y;
        Mesh[i].v2.z += Transform.z;
    }
}

void rotateMesh(std::vector<Tri4d_t>& Mesh, Vec4d_t Rotation) {//multiplies two matricies for rotation
    int i, j, k, l; //w doesnt change with rotation, no need to alter mult matrix since we dodge the 4th coordinate (w) with horrible pointer offset code
    Tri4d_t Rtri;
    float sinRZ = sin(Rotation.x);//alpha
    float cosRZ = cos(Rotation.x);
    float sinRY = sin(Rotation.y);//beta
    float cosRY = cos(Rotation.y);
    float sinRX = sin(Rotation.z);//gamma
    float cosRX = cos(Rotation.z);
    float multMatrix[3][3] = { {cosRZ * cosRY, cosRZ * sinRY * sinRX - sinRZ * cosRX                , cosRZ * sinRY * cosRX + sinRZ * sinRX}
                              ,{sinRZ * cosRY, sinRZ * sinRY * sinRX + cosRZ * cosRX        , sinRZ * sinRZ * sinRY * cosRX - cosRZ * sinRX}
                              ,{-sinRY       , cosRY * sinRX                                        , cosRY * cosRX                        }};

    for (i = 0; i < Mesh.size(); i++) {
        std::memcpy(&Rtri, &Mesh[i], sizeof(Tri4d_t));
        for (l = 0; l < 3; l++) {//iterate through triangle vertice pairs
            for (j = 0; j < 3; j++) {//dot product loop
                float sum = 0;
                for (k = 0; k < 3; k++) {//very safe memory accessing inbound
                    //so its accessing the pointer to the start of Mesh[] and adjusted the pointer offset by + ((...)) and reading the resulting address as a pointer to a float :)
                    sum += *(float*)((unsigned char*)&Mesh[i] + ((l * 16) + (k * 4))) * multMatrix[j][k];
                }//some more too
                //same thing here, fun fact, l is the y value of the matrix and k is the x of the matrix, i think ||| ITS NOT
                *(float*)((unsigned char*)&Rtri + ((l * 16) + (j * 4))) = sum;
            }
        }
        std::memcpy(&Mesh[i], &Rtri, sizeof(Tri4d_t));
    }
}

void worldSpaceMesh(std::vector<Tri4d_t>& Mesh, MeshInstance_t Instance) {//move a mesh based on its own position and rotation as an instance of a parent mesh
    rotateMesh(Mesh, Instance.rot);
    translateMesh(Mesh, Instance.pos);
    
}

void cameraSpaceMesh(std::vector<Tri4d_t>& Mesh, Camera_t camera) {//move the inverted position and rotation to the world scene
    Vec4d_t tPos = { -camera.position.x, -camera.position.y, -camera.position.z, -camera.position.w };
    Vec4d_t tRot = { -camera.rotation.x, -camera.rotation.y, -camera.rotation.z, -camera.rotation.w };
    translateMesh(Mesh, tPos);
    rotateMesh(Mesh, tRot);    
}

int clipOrCull(MeshInstance_t mesh, Camera_t camera) {//early rejection of meshes based on set bounding box coordinates and clip plane inequalities
  /*if fully in = 8
    if fully out = 0
    if parial = < 8 */
    int i;
    std::vector<Tri4d_t> tempBounds = mesh.parentMesh.boundBox;
    worldSpaceMesh(tempBounds, mesh);
    cameraSpaceMesh(tempBounds, camera);
    projectMeshMatrix(tempBounds);
    //we in clip space now
    int result = 0;
    for(i = 0; i < tempBounds.size(); i++){
        if (tempBounds[i].v0.w <= 0) result;//axe this soon
        if ((tempBounds[i].v0.x >= -tempBounds[i].v0.w) && (tempBounds[i].v0.x <= tempBounds[i].v0.w) &&
            (tempBounds[i].v0.y >= -tempBounds[i].v0.w) && (tempBounds[i].v0.y <= tempBounds[i].v0.w) &&
            (tempBounds[i].v0.z >= -tempBounds[i].v0.w) && (tempBounds[i].v0.z <= tempBounds[i].v0.w)) {
            result += 1;//above test for inside, inside on all planes -> point is inside
        }
    }
    return result;
}

void clipTriangle(std::vector<Tri4d_t>& output, Tri4d_t input, Plane_t plane) {//clip and redraw triangles along a clipping plane
    Vec4d_t A, B, C, TA, TB, TC;
    float d0, d1, d2;
    d0 = signedDistance(plane, input.v0);
    d1 = signedDistance(plane, input.v1);
    d2 = signedDistance(plane, input.v2);
    if (d0 >= 0 && d1 >= 0 && d2 >= 0) {//no clipping needed, do nothing            
        output.push_back(input);
    }
    else if (d0 < 0 && d1 < 0 && d2 < 0) {//all out of the selection      
        return;
    }
    else if (((d0 > 0) + (d1 > 0) + (d2 > 0)) == 1) {//only one positive
        //A will be the positive vertex
        if (d0 > d1) {
            if (d0 > d2) {
                A = input.v0;
                B = input.v1;
                C = input.v2;
            }
            else {
                A = input.v2;
                B = input.v1;
                C = input.v0;
            }
        }
        else {
            if (d1 > d2) {
                A = input.v1;
                B = input.v2;
                C = input.v0;
            }
            else {
                A = input.v2;
                B = input.v1;
                C = input.v0;
            }
        }
        TB = linePlaneIntersection(A, B, plane);
        TC = linePlaneIntersection(A, C, plane);
        output.push_back(Tri4d_t{ A, TB, TC });
    }
    else if (((d0 < 0) + (d1 < 0) + (d2 < 0)) == 1) {//only one negative
        //c will be the negative one
        if (d0 < d1) {
            if (d0 < d2) {
                C = input.v0;
                B = input.v1;
                A = input.v2;
            }
            else {
                C = input.v2;
                B = input.v1;
                A = input.v0;
            }
        }
        else {
            if (d1 < d2) {
                C = input.v1;
                B = input.v2;
                A = input.v0;
            }
            else {
                C = input.v2;
                B = input.v1;
                A = input.v0;
            }
        }
        TA = linePlaneIntersection(A, C, plane);
        TB = linePlaneIntersection(B, C, plane);
        output.push_back(Tri4d_t{ A, B, TA });
        output.push_back(Tri4d_t{ TA, B, TB });
    }
}

std::vector<Tri4d_t> clipMesh(std::vector<Tri4d_t>& mesh, Camera_t camera) {//handle passing a mesh through the clipper along all its planes without any data hazards
    int i, j;
    Vec4d_t A, B, C, TA, TB, TC;
    std::vector<Tri4d_t> cInput = mesh;
    std::vector<Tri4d_t> cOutput;
    for (j = 0; j < camera.clippingPlanes.size(); j++) {
        cOutput.clear();
        for (i = 0; i < cInput.size(); i++) {
            clipTriangle(cOutput, cInput[i], camera.clippingPlanes[j]);
        }
        cInput = cOutput;
    }
    return cInput;
}

std::vector<Tri4d_t> rotTranClipMesh(struct MeshInstance_t mesh, struct Camera_t camera ) {//handles the entire clipping pipeline with a mesh and a camera
    int i, cResult;
    std::vector<Tri4d_t> clipped;
    cResult = clipOrCull(mesh, camera);
    if (cResult == 0) {//fully out
        std::vector<Tri4d_t> empty;
        return empty;
    }
    std::vector<Tri4d_t> processing = mesh.parentMesh.Triangles;
    worldSpaceMesh(processing, mesh);
    cameraSpaceMesh(processing, camera);
    projectMeshMatrix(processing);
    if (cResult < 8) {//partial
        clipped = clipMesh(processing, camera);
    }
    else {//fully in
        clipped = processing;
    }
    return clipped;
}

void drawLine(float x0, float y0, float x1, float y1) {//bresenham line drawing between two points
    float x, y, dx, dy, step;
    int i;
    dx = (x1 - x0);
    dy = (y1 - y0);
    if (fabs(dx) >= fabs(dy)) {
        step = fabs(dx);
    }
    else {
        step = fabs(dy);
    }
    dx = dx / step;
    dy = dy / step;
    x = x0;
    y = y0;
    i = 0;
    int xmod, ymod;
    while (i <= step) {
        int xr = round(x);
        int yr = round(y);  
        setPixel(xr, yr, 0xffffffff);
        x = x + dx;
        y = y + dy;
        i++;
    }
}

void drawWireTri(Tri2d_t Tri) {//handles drawing three lines between verticies of a triangle
    drawLine(Tri.v0.x, Tri.v0.y, Tri.v1.x, Tri.v1.y);
    drawLine(Tri.v1.x, Tri.v1.y, Tri.v2.x, Tri.v2.y);
    drawLine(Tri.v2.x, Tri.v2.y, Tri.v0.x, Tri.v0.y);
}

void drawFullTri(Tri2d_t Tri, Uint32 Colour) {//draws a filled in triangle of a set colour, no zbuffering im lazy
    //swap so we get a nice y ordered tringle
    if (Tri.v1.y < Tri.v0.y) { swapVec2d_t(&Tri.v1, &Tri.v0); }
    if (Tri.v2.y < Tri.v0.y) { swapVec2d_t(&Tri.v2, &Tri.v0); }
    if (Tri.v2.y < Tri.v1.y) { swapVec2d_t(&Tri.v2, &Tri.v1); }
    std::vector<int> x01, x12, x02, x012, x_left, x_right;
    //generates three arrays of every x value that lies on the triangle edges
    interpolatei(x01, Tri.v0.y, Tri.v0.x, Tri.v1.y, Tri.v1.x);
    interpolatei(x12, Tri.v1.y, Tri.v1.x, Tri.v2.y, Tri.v2.x);
    interpolatei(x02, Tri.v0.y, Tri.v0.x, Tri.v2.y, Tri.v2.x);
    x01.pop_back();
    x012.reserve(x01.size() + x12.size());
    x012.insert(x012.end(), x01.begin(), x01.end());
    x012.insert(x012.end(), x12.begin(), x12.end());
    //find which is to the left or right
    int m = floor(x02.size() / 2);
    if (x02[m] < x012[m]) {
        x_left = x02;
        x_right = x012;
    }
    else {
        x_left = x012;
        x_right = x02;
    }
    //actual drawing
    int y, x, xmod, ymod;
    for (y = Tri.v0.y; y < Tri.v2.y; y++) {
        for (x = x_left[y - Tri.v0.y]; x < x_right[y - Tri.v0.y]; x++) {
            setPixel(x, y, Colour);
        }
    }
    x01.clear();
    x02.clear();
    x12.clear();
    x012.clear();
    x_left.clear();
    x_right.clear();
}
 
void drawShadedTri(Tri2d_t Tri, int Colour) {//draws a triangle with interpolated colour values corresponding to vx.h
    //swap so we get a nice y ordered tringle
    if (Tri.v1.y < Tri.v0.y) { swapVec2d_t(&Tri.v1, &Tri.v0); }
    if (Tri.v2.y < Tri.v0.y) { swapVec2d_t(&Tri.v2, &Tri.v0); }
    if (Tri.v2.y < Tri.v1.y) { swapVec2d_t(&Tri.v2, &Tri.v1); }
    std::vector<int> x01, x12, x02, x012, x_left, x_right;
    std::vector<float> h01, h12, h02, h012, h_segment, h_right, h_left;
    std::vector<double> z01, z12, z02, z012, z_segment, z_right, z_left;
    //generates three arrays of every x value that lies on the triangle edges
    interpolatei(x01, Tri.v0.y, Tri.v0.x, Tri.v1.y, Tri.v1.x);
    interpolatef(h01, Tri.v0.y, Tri.v0.h, Tri.v1.y, Tri.v1.h);
    interpolated(z01, Tri.v0.y, Tri.v0.z, Tri.v1.y, Tri.v1.z);

    interpolatei(x12, Tri.v1.y, Tri.v1.x, Tri.v2.y, Tri.v2.x);
    interpolatef(h12, Tri.v1.y, Tri.v1.h, Tri.v2.y, Tri.v2.h);
    interpolated(z12, Tri.v1.y, Tri.v1.z, Tri.v2.y, Tri.v2.z);

    interpolatei(x02, Tri.v0.y, Tri.v0.x, Tri.v2.y, Tri.v2.x);
    interpolatef(h02, Tri.v0.y, Tri.v0.h, Tri.v2.y, Tri.v2.h);
    interpolated(z02, Tri.v0.y, Tri.v0.z, Tri.v2.y, Tri.v2.z);

    if(x01.size() != 0){ x01.pop_back(); }
    if(h01.size() != 0){ h01.pop_back(); }
    if(z01.size() != 0){ z01.pop_back(); }
    
    x012.reserve(x01.size() + x12.size());
    x012.insert(x012.end(), x01.begin(), x01.end());
    x012.insert(x012.end(), x12.begin(), x12.end());

    h012.reserve(h01.size() + h12.size());
    h012.insert(h012.end(), h01.begin(), h01.end());
    h012.insert(h012.end(), h12.begin(), h12.end());

    z012.reserve(z01.size() + z12.size());
    z012.insert(z012.end(), z01.begin(), z01.end());
    z012.insert(z012.end(), z12.begin(), z12.end());

    //find which is to the left or right
    int m = floor(x02.size() / 2);
    if (x02[m] < x012[m]) {
        x_left = x02;
        h_left = h02;
        z_left = z02;
        x_right = x012;
        h_right = h012;
        z_right = z012;
    }
    else {
        x_left = x012;
        h_left = h012;
        z_left = z012;
        x_right = x02;
        h_right = h02;
        z_right = z02;
    }
    //actual drawing
    int y, x, xmod, ymod, x_l, x_r, shaded_colour;
    for (y = Tri.v0.y; y < Tri.v2.y; y++) {
        x_l = x_left[y - Tri.v0.y];
        x_r = x_right[y - Tri.v0.y];
        interpolatef(h_segment, x_l, h_left[y - Tri.v0.y], x_r, h_right[y - Tri.v0.y]);
        interpolated(z_segment, x_l, z_left[y - Tri.v0.y], x_r, z_right[y - Tri.v0.y]);
        y = y;
        for (x = x_l; x <= x_r; x++) {
            shaded_colour = multiplyColour(Colour, h_segment[x - x_l]);
            shaded_colour = shaded_colour | 0x000000ff;
            setZPixel(x, y, z_segment[x - x_l], shaded_colour);
        }
        h_segment.clear();
    }
    x01.clear(); x12.clear(); x02.clear(); x012.clear(); x_left.clear(); x_right.clear();
    h02.clear(); h12.clear(); h02.clear(); h012.clear(); h_left.clear(); h_right.clear();
    z02.clear(); z12.clear(); z02.clear(); z012.clear(); z_left.clear(); z_right.clear();
}

bool update()
{
    SDL_Event e;
    if (SDL_PollEvent(&e))
    {
        if (e.type == SDL_EVENT_QUIT)
        {
            return false;
        }
        if (e.type == SDL_EVENT_KEY_UP && e.key.key == SDLK_ESCAPE)
        {
            return false;
        }
    }

    char* pix;
    int pitch;

    SDL_LockTexture(gSDLTexture, NULL, (void**)&pix, &pitch);
    for (int i = 0, sp = 0, dp = 0; i < WINDOW_HEIGHT; i++, dp += WINDOW_WIDTH, sp += pitch)
        memcpy(pix + sp, gFrameBuffer + dp, WINDOW_WIDTH * 4);

    SDL_UnlockTexture(gSDLTexture);
    SDL_RenderTexture(gSDLRenderer, gSDLTexture, NULL, NULL);
    SDL_RenderPresent(gSDLRenderer);
    SDL_Delay(1);//this eats 3 frames per second : (
    return true;
}

float Loop = 1;

void static setup_scene() {//ran once to mesh meshes and instance instances

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

    //struct MeshInstance_t mesh1Instance3;
    //mesh1Instance3.parentMesh = cubeCube;
    //mesh1Instance3.pos = { 150,0,0 };
    

    struct Camera_t camera;
    camera.position = { 0,30,300 };//by default facing down -Z
    camera.rotation = { 0,0,0 };

    rendererScene.meshInstances.push_back(mesh1Instance1);
    rendererScene.meshInstances.push_back(mesh1Instance2);
    //rendererScene.meshInstances.push_back(mesh1Instance3);
    rendererScene.meshInstanceCount = rendererScene.meshInstances.size();
    rendererScene.camera = camera;

}

void render_scene(struct Scene_t scene) {//rendering pipeline brought together
    int i, j;
    std::vector<Tri3d_t> tempProjected;
    for (i = 0; i < scene.meshInstanceCount; i++) {
        std::vector<Tri4d_t> tempTri = rotTranClipMesh(scene.meshInstances[i], scene.camera);
        tempProjected = perspectiveDivide(tempTri);
        std::vector<Tri2d_t> tempScreen = viewportTransformation(tempProjected);
        for (j = 0; j < tempScreen.size(); j++) {
            drawShadedTri(tempScreen[j], (multiplyColour(0x110000ff, j)));
        }
        for (j = 0; j < tempScreen.size(); j++) {
            //drawWireTri(tempScreen[j]);
        }
        tempTri.clear(); tempProjected.clear(); tempScreen.clear();
    }
}


void render(Uint64 aTicks)/*does the funny rendering*/
{
    for (int i = 0, c = 0; i < WINDOW_HEIGHT; i++)//    bg black and zbuffer
    {
        for (int j = 0; j < WINDOW_WIDTH; j++, c++)
        {
            gFrameBuffer[c] = 0x000000ff;
            zBuffer[c] = -FAR;
        }
    }
    rendererScene.meshInstances[0].rot.x += .01;
    rendererScene.meshInstances[0].rot.y += .01;
    rendererScene.meshInstances[0].rot.z += .01;
    rendererScene.meshInstances[1].rot.y += .01;
    //rendererScene.camera.rotation.y += 0.01;
    render_scene(rendererScene);
   
}

void loop()
{
    if (!update())
    {
        gDone = 1;
#ifdef __EMSCRIPTEN__
        emscripten_cancel_main_loop();
#endif
    }
    else
{
        render(SDL_GetTicks());
    }
}

int main(int argc, char** argv)
{
    if (!SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS))
    {
        return -1;
    }

    gFrameBuffer = new int[WINDOW_WIDTH * WINDOW_HEIGHT];
    zBuffer = new double[WINDOW_WIDTH * WINDOW_HEIGHT];
    gSDLWindow = SDL_CreateWindow("SDL3 window", WINDOW_WIDTH, WINDOW_HEIGHT, 0);
    gSDLRenderer = SDL_CreateRenderer(gSDLWindow, NULL);
    gSDLTexture = SDL_CreateTexture(gSDLRenderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, WINDOW_WIDTH, WINDOW_HEIGHT);

    if (!gFrameBuffer || !gSDLWindow || !gSDLRenderer || !gSDLTexture)
        return -1;

    gDone = 0;
#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(loop, 0, 1);
#else
    setup_scene();
    while (!gDone)
    {
        loop();
    }
#endif
    SDL_DestroyTexture(gSDLTexture);
    SDL_DestroyRenderer(gSDLRenderer);
    SDL_DestroyWindow(gSDLWindow);
    SDL_Quit();

    return 0;
}