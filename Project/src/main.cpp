#include "../include/Grid.hpp"
// #include "rectangular.hpp"
#include "../include/Particle.hpp"
#include "../include/sphere.hpp"
#include "../include/triangle.hpp"
#include "../include/camera.hpp"
//#include "extra.hpp"
#include "../include/simDomain.hpp"
#include <GL/glut.h>
#include <iostream>

namespace
{
// Global variables here.

// This is the camera
Camera camera;

// These are state variables for the UI
bool gMousePressed = false;

// These are arrays for display lists
GLuint gAxisList;
GLuint gPointList;
GLuint gSurfaceLists;
GLuint gBBoxLineLists;

// One big snow particle group
SnowParticleSet* globalSPS;
Grid* globalGridMesh;
SimDomain* globalSimDomain;
bool printTestMesh = false; 
/* Do not have methods to load obj file now */
//MeshTriangle aTestMesh("../media/spot_triangulated_good.obj");

// Declarations of functions whose implementations occur later.
void arcballRotation(int endX, int endY);
void keyboardFunc(unsigned char key, int x, int y);
// void specialFunc(int key, int x, int y);
void mouseFunc(int button, int state, int x, int y);
void motionFunc(int x, int y);
void reshapeFunc(int w, int h);
void drawScene(void);
void initRendering();
// void loadObjects(int argc, char* argv[]);
void makeDisplayLists();
void updateParticleLists();

// I dont know how to run sim while in glut main loop
// trying
// bool runningNow = false;
float globalSimTime = 0;
void simulation();
bool started = false;

// This function is called whenever a "Normal" key press is
// received.
void keyboardFunc(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27:  // Escape key
            exit(0);
            break;
        case ' ':
        {
            Matrix4f eye = Matrix4f::Identity();
            camera.SetRotation(eye);
            camera.SetCenter(Vector3f(0, 0, 0));
            break;
        }
        case 's':
        case 'S':
            std::cout << "start simulation" << std::endl;
            // runningNow = true;
            simulation();
            break;
        case 'h':
            printTestMesh = !printTestMesh;
            break;
        default:
            std::cout << "Unhandled key press " << key << "." << std::endl;
    }

    glutPostRedisplay();
}

//  Called when mouse button is pressed.
void mouseFunc(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN)
    {
        gMousePressed = true;

        switch (button)
        {
            case GLUT_LEFT_BUTTON:
                camera.MouseClick(Camera::LEFT, x, y);
                break;
            case GLUT_MIDDLE_BUTTON:
                camera.MouseClick(Camera::MIDDLE, x, y);
                break;
            case GLUT_RIGHT_BUTTON:
                camera.MouseClick(Camera::RIGHT, x, y);
            default:
                break;
        }
    }
    else
    {
        camera.MouseRelease(x, y);
        gMousePressed = false;
    }
    glutPostRedisplay();
}

// Called when mouse is moved while button pressed.
void motionFunc(int x, int y)
{
    camera.MouseDrag(x, y);

    glutPostRedisplay();
}

// Called when the window is resized
// w, h - width and height of the window in pixels.
void reshapeFunc(int w, int h)
{
    camera.SetDimensions(w, h);

    camera.SetViewport(0, 0, w, h);
    camera.ApplyViewport();

    // Set up a perspective view, with square aspect ratio
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    camera.SetPerspective(50);
    camera.ApplyPerspective();
}

// This function is responsible for displaying the object.
void drawScene(void)
{
    // std::cout << "drawing " << std::endl;
    // Clear the rendering window
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Light color (RGBA)
    GLfloat Lt0diff[] = {1.0, 1.0, 1.0, 10.0};
    GLfloat Lt0pos[] = {0.0, 2.5, -2.5, 1.0};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, Lt0diff);
    glLightfv(GL_LIGHT0, GL_POSITION, Lt0pos);

    camera.ApplyModelview();

    // xyz axis
    glPushMatrix();
    glTranslated(camera.GetCenter()[0], camera.GetCenter()[1],
                 camera.GetCenter()[2]);
    glCallList(gAxisList);
    glPopMatrix();

    glCallList(gSurfaceLists);
    glCallList(gBBoxLineLists);
    glCallList(gPointList);

    // Dump the image to the screen.
    glutSwapBuffers();
}

// Initialize OpenGL's rendering modes
void initRendering()
{
    glEnable(GL_DEPTH_TEST);  // Depth testing must be turned on
    glEnable(GL_LIGHTING);    // Enable lighting calculations
    glEnable(GL_LIGHT0);      // Turn on light #0.

    // Setup polygon drawing
    glShadeModel(GL_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    // point size
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    // Clear to black
    glClearColor(52. / 255, 87. / 255, 115. / 255, 1);

    // Base material colors (they don't change)
    GLfloat diffColor[] = {0.4, 0.4, 0.4, 1};
    GLfloat specColor[] = {0.9, 0.9, 0.9, 1};
    GLfloat shininess[] = {50.0};

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, diffColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

void makeDisplayLists()
{
    // std::cout << " re-generating glut lists" << std::endl;
    gAxisList = glGenLists(1);
    gPointList = glGenLists(1);
    gSurfaceLists = glGenLists(1);
    gBBoxLineLists = glGenLists(1);

    // Compile the display lists

    glNewList(gSurfaceLists, GL_COMPILE);
    {
        // This will use the current material color and light
        // positions.  Just set these in drawScene();
        glEnable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // This tells openGL to *not* draw backwards-facing triangles.
        // This is more efficient, and in addition it will help you
        // make sure that your triangles are drawn in the right order.
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glBegin(GL_TRIANGLES);
        // make the floor
        // floor is a big tri-pair whose y is the min of global bbox
        float y = globalGridMesh->global_bbox.pMin.y();
        Vector3f normal(0, 1., 0);
        Vector3f fP0(-10., y, -10);
        Vector3f fP1(10., y, -10);
        Vector3f fP2(-10., y, 10);
        Vector3f fP3(10., y, 10);
        glNormal(normal);
        glVertex(fP0);
        glNormal(normal);
        glVertex(fP2);
        glNormal(normal);
        glVertex(fP1);
        glNormal(normal);
        glVertex(fP3);
        glNormal(normal);
        glVertex(fP1);
        glNormal(normal);
        glVertex(fP2);
        glEnd();
        glPopAttrib();
    }
    glEndList();

    glNewList(gAxisList, GL_COMPILE);
    {
        // Save current state of OpenGL
        glPushAttrib(GL_ALL_ATTRIB_BITS);

        // This is to draw the axes when the mouse button is down
        glDisable(GL_LIGHTING);
        glLineWidth(3);
        glPushMatrix();
        glScaled(5.0, 5.0, 5.0);
        glBegin(GL_LINES);
        glColor4f(1, 0.5, 0.5, 1);
        glVertex3d(0, 0, 0);
        glVertex3d(1, 0, 0);
        glColor4f(0.5, 1, 0.5, 1);
        glVertex3d(0, 0, 0);
        glVertex3d(0, 1, 0);
        glColor4f(0.5, 0.5, 1, 1);
        glVertex3d(0, 0, 0);
        glVertex3d(0, 0, 1);

        glColor4f(0.5, 0.5, 0.5, 1);
        glVertex3d(0, 0, 0);
        glVertex3d(-1, 0, 0);
        glVertex3d(0, 0, 0);
        glVertex3d(0, -1, 0);
        glVertex3d(0, 0, 0);
        glVertex3d(0, 0, -1);

        glEnd();
        glPopMatrix();

        glPopAttrib();
    }
    glEndList();

    glNewList(gPointList, GL_COMPILE);
    {
        // Save current state of OpenGL
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        // Setup for point drawing
        glDisable(GL_LIGHTING);

        float pointScale = 138. / globalSPS->particles[0]->m->lNumDensity;

        for (auto& oneSP : globalSPS->particles)
        {
            float relRho = oneSP->density / oneSP->m->initialDensity;
            float contrast = 0.6;
            relRho = relRho * contrast + 1. - contrast;
            float relVol = std::pow(
                oneSP->volume / globalGridMesh->eachCellVolume, 1. / 3.);
            // std::cout << " what is rel vol" << relVol << std::endl;
            // std::cout << " what is rel rho" << relRho << std::endl;
            // glColor4f(relRho, relRho, relRho, 1);
            glColor4f(relRho, relRho, relRho, 1);
            glPointSize(pointScale * relVol);
            glLineWidth(1);
            glBegin(GL_POINTS);
            glVertex(oneSP->position);
            glEnd();
        }

        glPopAttrib();
    }
    glEndList();

    glNewList(gBBoxLineLists, GL_COMPILE);
    {
        // Save current state of OpenGL
        glPushAttrib(GL_ALL_ATTRIB_BITS);

        // This is to draw the axes when the mouse button is down
        glDisable(GL_LIGHTING);
        glLineWidth(1);
        glPushMatrix();
        // glScaled(5.0, 5.0, 5.0);

        glBegin(GL_LINES);
        glColor4f(1, 1, 1, 0.5);
        Vector3f p0 = globalSimDomain->gridMesh->global_bbox.pMin;
        Vector3f p6 = globalSimDomain->gridMesh->global_bbox.pMax;
        Vector3f p1 = p0, p2 = p0, p3 = p0;
        p1.x() = p6.x();
        p3.y() = p6.y();
        p2.x() = p6.x();
        p2.y() = p6.y();
        Vector3f p4 = p6, p5 = p6, p7 = p6;
        p5.y() = p0.y();
        p7.x() = p0.x();
        p4.x() = p0.x();
        p4.y() = p0.y();

        glVertex(p0);
        glVertex(p1);
        glVertex(p1);
        glVertex(p2);
        glVertex(p2);
        glVertex(p3);
        glVertex(p3);
        glVertex(p0);

        glVertex(p4);
        glVertex(p5);
        glVertex(p5);
        glVertex(p6);
        glVertex(p6);
        glVertex(p7);
        glVertex(p7);
        glVertex(p4);

        glVertex(p0);
        glVertex(p4);
        glVertex(p1);
        glVertex(p5);
        glVertex(p2);
        glVertex(p6);
        glVertex(p3);
        glVertex(p7);

        glEnd();
        glPopMatrix();

        glPopAttrib();
    }
    glEndList();
}

void updateParticleLists()
{
    gPointList = glGenLists(1);
    glNewList(gPointList, GL_COMPILE);
    {
        // Save current state of OpenGL
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        // Setup for point drawing
        glDisable(GL_LIGHTING);

        float pointScale = 138. / globalSPS->particles[0]->m->lNumDensity;

        for (auto& oneSP : globalSPS->particles)
        {
            float relRho = oneSP->density / oneSP->m->initialDensity;
            float contrast = 0.6;
            relRho = relRho * contrast + 1. - contrast;
            float relVol = std::pow(
                oneSP->volume / globalGridMesh->eachCellVolume, 1. / 3.);
            // std::cout << " what is rel vol" << relVol << std::endl;
            // std::cout << " what is rel rho" << relRho << std::endl;
            // glColor4f(relRho, relRho, relRho, 1);
            glColor4f(relRho, relRho, relRho, 1);
            glPointSize(pointScale * relVol);
            glLineWidth(1);
            glBegin(GL_POINTS);
            glVertex(oneSP->position);
            glEnd();
        }

        glPopAttrib();
    }
    glEndList();
}

void simulation()
{
    // init
    float endTime = globalSimTime + 1. / 5.;
    while (globalSimTime < endTime)
    {
        std::cout << " running now, time is " << globalSimTime << std::endl;
        globalSimTime += deltaT;
        // test
        // globalSPS->update();
        // run real simulation onetime here
        globalSimDomain->oneTimeSimulate();
        if ((deltaT >= 1. / FRAMERATE) ||
            (std::abs(std::remainder(globalSimTime, 1. / FRAMERATE)) <
             0.5 * deltaT))
        {
            // should re-draw now
            std::cout << "re-drawing" << std::endl;
            updateParticleLists();
            drawScene();
        }
    }
    std::cout << " simulation paused at " << globalSimTime << std::endl;
    return;
}

}  // namespace

// Main function for snow simulation
int main(int argc, char** argv)
{
     // snow sim
    SnowParticleMaterial m;
    m.lNumDensity = 35;
    //MeshTriangle cow("../media/spot_triangulated_good.obj");
    globalSPS = new SnowParticleSet();
    //globalSPS->addParticlesInAShape(&cow, 10. * Vector3f(0.2, 0, -1.), &m);
    //SnowParticleSet mirroredCow;
    //mirroredCow.CreateMirror(*globalSPS, 0, 0, 1., 2.5, Vector3f(0, 0, -2.5));
    //globalSPS->appendSet(mirroredCow);
    Bounds3 bbox(Vector3f(-1., 10, -7.), Vector3f(1., 0, 1.5));
    //bbox = Union(bbox, cow.getBounds());
    // add a carpet of snow using the bbox and rectangle
    //Vector3f floorP0(bbox.pMin);
    //Vector3f floorP1(bbox.pMax);
    //floorP1.y() = floorP0.y() - 0.06;
    //Rectangular floor(floorP0, floorP1);
    //globalSPS->addParticlesInAShape(&floor, &m);
    Vector3f center = (bbox.pMin + bbox.pMax)/3.0;
    float radius = 0.5;
    Sphere sphere(center,radius);
    globalSPS->addParticlesInAShape(&sphere, Vector3f(0, -50, 0), &m);
    //bbox = Union(bbox, floor.getBounds());
    bbox = Union(bbox, sphere.getBounds());

    // Mesh grid and simulation domain
    globalGridMesh = new Grid(bbox, Vector3f(.05, .05, .05), globalSPS);
    globalSimDomain = new SimDomain(globalSPS, globalGridMesh);
    globalSimDomain->initializeSimulator();

    glutInit(&argc, argv);
    // We're going to animate it, so double buffer
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    // Initial parameters for window position and size
    glutInitWindowPosition(60, 60);
    ;
    glutInitWindowSize(600, 600);
    // set initial camera
    camera.SetDimensions(600, 600);
    camera.SetDistance(10);
    camera.SetCenter(Vector3f(0, 0, 0));
    // create window
    glutCreateWindow("Snow Simulation");
    // Initialize OpenGL parameters.
    initRendering();
    // Set up callback functions for key presses
    glutKeyboardFunc(keyboardFunc);
    // Set up callback functions for mouse
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);
    // Set up the callback function for resizing windows
    glutReshapeFunc(reshapeFunc);
    // Call this whenever window needs redrawing
    glutDisplayFunc(drawScene);
    // compiled gl list can be drawed at convenience later
    makeDisplayLists();
    // Start the main loop.  glutMainLoop never returns.
    glutMainLoop();

    delete (globalSimDomain);
    delete (globalGridMesh);
    delete (globalSPS);
    return 0;
}