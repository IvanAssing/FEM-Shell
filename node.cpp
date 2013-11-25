#include "node.h"
#include <QGLWidget>


Node::Node()
{
}


Node::Node(int index_, double x_, double y_)
    :index(index_), x(x_), y(y_)
{

}


void Node::draw(void)
{
    glColor4d(1.0, 0.0, 0.0, 0.5);
    glPointSize(5.f);
    glBegin(GL_POINTS);
    glVertex3d(x, y, 0.0);
    glEnd();
}
