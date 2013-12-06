#include "node.h"
#include <QGLWidget>

#include <cmath>


Node::Node()
{
}


Node::Node(int index_, double x_, double y_, double z_)
    :index(index_), x(x_), y(y_), z(z_)
{
    loadValues = new double[6];
    lockStatus = new bool[6];

    for(int i=0; i<6; i++)
    {
        loadValues[i] = 0.0;
        lockStatus[i] = false;
    }

}

void Node::setup(Boundary &b)
{
    for(int i=0; i<6; i++)
    {
        loadValues[i] = b.loadValues[i];
        lockStatus[i] = b.lockStatus[i];
    }
}


void Node::draw(void)
{
    glColor4d(1.0, 1.0, 1.0, 1.0);

    glPointSize(8.f);
    glBegin(GL_POINTS);
    glVertex3d(x, y, z);
    glEnd();
}

void Node::draw_lock(void)
{
    glColor4d(1.0, 1.0, 1.0, 0.8);
    glLineWidth(5.0f);

    double dt = 0.01;
    if(lockStatus[0])
    {
        glBegin(GL_LINES);
        glVertex3d(x-dt, y, z);
        glVertex3d(x+dt, y, z);
        glEnd();
    }
    if(lockStatus[1])
    {
        glBegin(GL_LINES);
        glVertex3d(x, y-dt, z);
        glVertex3d(x, y+dt, z);
        glEnd();
    }
    if(lockStatus[2])
    {
        glBegin(GL_LINES);
        glVertex3d(x, y, z-dt);
        glVertex3d(x, y, z+dt);
        glEnd();
    }

    double ds = dt/2.;
    int n = 50;
    double ns = 2*M_PI/n;
    if(lockStatus[3])
    {
        glBegin(GL_LINE_LOOP);
        for(int i=0; i<n; i++)
            glVertex3d(x, y + ds*cos(i*ns), z + ds*sin(i*ns));
        glEnd();
    }
    if(lockStatus[4])
    {
        glBegin(GL_LINE_LOOP);
        for(int i=0; i<n; i++)
            glVertex3d(x + ds*cos(i*ns), y, z + ds*sin(i*ns));
        glEnd();
    }
    if(lockStatus[5])
    {
        glBegin(GL_LINE_LOOP);
        for(int i=0; i<n; i++)
            glVertex3d(x + ds*cos(i*ns), y + ds*sin(i*ns), z);
        glEnd();
    }

}

void Node::draw_load(void)
{

    glPointSize(10.0f);
    glColor4d(1.0, 0.0, 0.8, 1.0);

    for(int i=0; i<6; i++)
        if(fabs(loadValues[i])> 1.0e-5)
{
            glBegin(GL_POINTS);
            glVertex3d(x, y, z);
            glEnd();
        }


}

Node::~Node()
{

}
