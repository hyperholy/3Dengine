#include "graphics_structs.h"
#include "render.h"
#include "scene.h"
#include "SDL3/SDL_main.h"

static int gDone;
SDL_Window* gSDLWindow;
SDL_Renderer* gSDLRenderer; 
SDL_Texture* gSDLTexture;

/*global variable scene*/
Scene_t rendererScene;

/*
TODO
30fps -> 90fps
use bitshifts for colour mult
better early rejection for h_segment where its all 1.0
occulision culling
shading
textures
*/

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
    SDL_RenderClear(gSDLRenderer);
    SDL_Delay(1);//this eats 7 frames per second : (
    return true;
}

void render(Uint64 aTicks)/*does the funny rendering*/
{
    rendererScene.meshInstances[0].rot.x += .01;
    rendererScene.meshInstances[0].rot.y += .01;
    rendererScene.meshInstances[0].rot.z += .01;
    //rendererScene.meshInstances[1].rot.y += .01;
    rendererScene.camera.rotation.y += 0.01;
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
    gSDLWindow = SDL_CreateWindow("SDL3 window", WINDOW_WIDTH, WINDOW_HEIGHT, SDL_WINDOW_FULLSCREEN);
    gSDLRenderer = SDL_CreateRenderer(gSDLWindow, NULL);
    gSDLTexture = SDL_CreateTexture(gSDLRenderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, WINDOW_WIDTH, WINDOW_HEIGHT);

    if (!gFrameBuffer || !gSDLWindow || !gSDLRenderer || !gSDLTexture)
        return -1;

    gDone = 0;
#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(loop, 0, 1);
#else
    setup_scene(rendererScene);
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