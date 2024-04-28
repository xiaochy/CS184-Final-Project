#include <cmath>
// #include <glad/glad.h>

#include <CGL/vector3D.h>
#include <nanogui/nanogui.h>

#include "../include/simulator.h"

#include "camera.h"
#include "camera_info.h"
// // Needed to generate stb_image binaries. Should only define in exactly one source file importing stb_image.h.
// #define STB_IMAGE_IMPLEMENTATION
// #include "misc/stb_image.h"

using namespace nanogui;
using namespace std;

void load(Grid& grid) {
    this->grid = grid;
}

void drawContents() {

}

// ----------------------------------------------------------------------------
// EVENT HANDLING
// ----------------------------------------------------------------------------

bool SnowSimulator::cursorPosCallbackEvent(double x, double y)
{
    if (left_down && !middle_down && !right_down)
    {
        if (ctrl_down)
        {
            mouseRightDragged(x, y);
        }
        else
        {
            mouseLeftDragged(x, y);
        }
    }
    else if (!left_down && !middle_down && right_down)
    {
        mouseRightDragged(x, y);
    }
    else if (!left_down && !middle_down && !right_down)
    {
        mouseMoved(x, y);
    }

    mouse_x = x;
    mouse_y = y;

    return true;
}

bool SnowSimulator::mouseButtonCallbackEvent(int button, int action,
                                             int modifiers)
{
    switch (action)
    {
    case GLFW_PRESS:
        switch (button)
        {
        case GLFW_MOUSE_BUTTON_LEFT:
            left_down = true;
            break;
        case GLFW_MOUSE_BUTTON_MIDDLE:
            middle_down = true;
            break;
        case GLFW_MOUSE_BUTTON_RIGHT:
            right_down = true;
            break;
        }
        return true;

    case GLFW_RELEASE:
        switch (button)
        {
        case GLFW_MOUSE_BUTTON_LEFT:
            left_down = false;
            break;
        case GLFW_MOUSE_BUTTON_MIDDLE:
            middle_down = false;
            break;
        case GLFW_MOUSE_BUTTON_RIGHT:
            right_down = false;
            break;
        }
        return true;
    }

    return false;
}

void SnowSimulator::mouseMoved(double x, double y) { y = screen_h - y; }

void SnowSimulator::mouseLeftDragged(double x, double y)
{
    float dx = x - mouse_x;
    float dy = y - mouse_y;

    camera.rotate_by(-dy * (PI / screen_h), -dx * (PI / screen_w));
}

void SnowSimulator::mouseRightDragged(double x, double y)
{
    camera.move_by(mouse_x - x, y - mouse_y, canonical_view_distance);
}

bool SnowSimulator::keyCallbackEvent(int key, int scancode, int action,
                                     int mods)
{
    ctrl_down = (bool)(mods & GLFW_MOD_CONTROL);

    if (action == GLFW_PRESS)
    {
        switch (key)
        {
        case GLFW_KEY_ESCAPE:
            is_alive = false;
            break;
        // TODO: we don't support reset operations rn
        // case 'r':
        // case 'R':
        //     grid->reset();
        //     break;
        case ' ':
            resetCamera();
            break;
        case 'p':
        case 'P':
            is_paused = !is_paused;
            break;
        case 'n':
        case 'N':
            if (is_paused)
            {
                is_paused = false;
                drawContents();
                is_paused = true;
            }
            break;
        }
    }

    return true;
}

bool SnowSimulator::dropCallbackEvent(int count, const char **filenames)
{
    return true;
}

bool SnowSimulator::scrollCallbackEvent(double x, double y)
{
    camera.move_forward(y * scroll_rate);
    return true;
}

bool SnowSimulator::resizeCallbackEvent(int width, int height)
{
    screen_w = width;
    screen_h = height;

    camera.set_screen_size(screen_w, screen_h);
    return true;
}

void SnowSimulator::mouseMoved(double x, double y) { y = screen_h - y; }

void SnowSimulator::mouseLeftDragged(double x, double y) {
  float dx = x - mouse_x;
  float dy = y - mouse_y;

  camera.rotate_by(-dy * (PI / screen_h), -dx * (PI / screen_w));
}

void SnowSimulator::mouseRightDragged(double x, double y) {
  camera.move_by(mouse_x - x, y - mouse_y, canonical_view_distance);
}