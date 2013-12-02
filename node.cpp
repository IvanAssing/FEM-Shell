#include "node.h"
#include <QGLWidget>


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
    glPointSize(3.f);
    glBegin(GL_POINTS);
    glVertex3d(x, y, 0.0);
    glEnd();
}

Node::~Node()
{

}
