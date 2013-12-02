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
    glColor4d(1.0, 1.0, 1.0, 0.9);

    glPointSize(3.f);
    glBegin(GL_POINTS);
    glVertex3d(x, y, z);
    glEnd();
}

void Node::draw_lock(void)
{
    glColor4d(1.0, 1.0, 1.0, 0.8);
    glLineWidth(5.0f);

    double dt = 0.1;
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
    int n = 20;
    double ns = 2*M_PI/20;
    if(lockStatus[3])
    {
        glBegin(GL_LINE_STRIP);
        for(int i=0; i<n; i++)
            glVertex3d(x, y + ds*cos(i*ns), z + ds*sin(i*ns));
        glEnd();
    }
    if(lockStatus[4])
    {
        glBegin(GL_LINE_STRIP);
        for(int i=0; i<n; i++)
            glVertex3d(x + ds*cos(i*ns), y, z + ds*sin(i*ns));
        glEnd();
    }
    if(lockStatus[5])
    {
        glBegin(GL_LINE_STRIP);
        for(int i=0; i<n; i++)
            glVertex3d(x + ds*cos(i*ns), y + ds*sin(i*ns), z);
        glEnd();
    }

}

void Node::draw_load(void)
{

    double dt = 0.1;
    double db = 0.8*dt;
    double dh = dt;
    double sg;
    dt *= 4;
    double tol = 1.0e-5;
    glColor4d(1.0, 0.0, 1.0, 0.8);
    glLineWidth(5.0f);

    if(fabs(loadValues[0])> tol)
   {
        sg = loadValues[0]>0.0 ? +1.0 : -1.0;
        glBegin(GL_LINES);
        glVertex3d(x, y, z);
        glVertex3d(x+sg*dt, y, z);
        glVertex3d(x, y, z);
        glVertex3d(x+sg*dh, y, z+db);
        glVertex3d(x, y, z);
        glVertex3d(x+sg*dh, y, z-db);
        glEnd();
    }
    if(fabs(loadValues[1])> tol)
   {
        sg = loadValues[1]>0.0 ? +1.0 : -1.0;
        glBegin(GL_LINES);
        glVertex3d(x, y, z);
        glVertex3d(x, y+sg*dt, z);
        glVertex3d(x, y, z);
        glVertex3d(x, y+sg*dh, z+db);
        glVertex3d(x, y, z);
        glVertex3d(x, y+sg*dh, z-db);
        glEnd();
    }
    if(fabs(loadValues[2])> tol)
   {
        sg = loadValues[2]>0.0 ? +1.0 : -1.0;
        glBegin(GL_LINES);
        glVertex3d(x, y, z);
        glVertex3d(x, y+sg*dt, z);
        glVertex3d(x, y, z);
        glVertex3d(x+db, y, z+sg*dh);
        glVertex3d(x, y, z);
        glVertex3d(x-db, y, z+sg*dh);
        glEnd();
    }

    if(fabs(loadValues[3])> tol)
   {
        sg = loadValues[3]>0.0 ? +1.0 : -1.0;
        glBegin(GL_LINES);
        glVertex3d(x, y, z);
        glVertex3d(x+sg*dt, y, z);
        glVertex3d(x, y, z);
        glVertex3d(x+sg*dh, y, z+db);
        glVertex3d(x, y, z);
        glVertex3d(x+sg*dh, y, z-db);
        glVertex3d(x+0.5*sg*dh, y, z);
        glVertex3d(x+1.5*sg*dh, y, z+db);
        glVertex3d(x+0.5*sg*dh, y, z);
        glVertex3d(x+1.5*sg*dh, y, z-db);
        glEnd();
    }
    if(fabs(loadValues[4])> tol)
   {
        sg = loadValues[4]>0.0 ? +1.0 : -1.0;
        glBegin(GL_LINES);
        glVertex3d(x, y, z);
        glVertex3d(x, y+sg*dt, z);
        glVertex3d(x, y, z);
        glVertex3d(x, y+sg*dh, z+db);
        glVertex3d(x, y, z);
        glVertex3d(x, y+sg*dh, z-db);
        glVertex3d(x, y+0.5*sg*dh, z);
        glVertex3d(x, y+1.5*sg*dh, z+db);
        glVertex3d(x, y+0.5*sg*dh, z);
        glVertex3d(x, y+1.5*sg*dh, z-db);
        glEnd();
    }
    if(fabs(loadValues[5])> tol)
   {
        sg = loadValues[5]>0.0 ? +1.0 : -1.0;
        glBegin(GL_LINES);
        glVertex3d(x, y, z);
        glVertex3d(x, y+sg*dt, z);
        glVertex3d(x, y, z);
        glVertex3d(x+db, y, z+sg*dh);
        glVertex3d(x, y, z);
        glVertex3d(x-db, y, z+sg*dh);
        glVertex3d(x, y, z+0.5*sg*dh);
        glVertex3d(x+db, y, z+1.5*sg*dh);
        glVertex3d(x, y, z+0.5*sg*dh);
        glVertex3d(x-db, y, z+1.5*sg*dh);
        glEnd();
    }

}

Node::~Node()
{

}
