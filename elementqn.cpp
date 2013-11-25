#include "elementqn.h"
#include <QGLWidget>
#include "matrix.h"
#include "polynomial2d.h"


ElementQN::ElementQN(Node **nodes_)
{
    nodes = new Node*[16];

    for(int i=0; i<16; i++)
        nodes[i] = nodes_[i];
}

void ElementQN::draw(void)
{

    glColor4d(0.0, 1.0, 0.0, 0.8);
    glLineWidth(2.5f);

    glBegin(GL_LINE_LOOP);
    {
        glVertex3d(nodes[0]->x, nodes[0]->y, 0.0);
        glVertex3d(nodes[3]->x, nodes[3]->y, 0.0);
        glVertex3d(nodes[15]->x, nodes[15]->y, 0.0);
        glVertex3d(nodes[12]->x, nodes[12]->y, 0.0);
    }
    glEnd();
}


void ElementQN::getStiffnessMatrix(Matrix &k, Polynomial2D **Bf, Polynomial2D **Bc)
{

//    for(int ii=0; ii<3; ii++)
//        for(int ij=0; ij<3; ij++)
//            for(int i=0; i<3; i++)
//                for(int j=0; j<3; j++)
//                    k.add(3*index[ii] + i, 3*index[ij] + j, (
//                          BtDB[3*ii+i][3*ij+j](0.5, 0.0) +
//                            BtDB[3*ii+i][3*ij+j](0.0, 0.5) +
//                            BtDB[3*ii+i][3*ij+j](0.5, 0.5))*_2A_by3);
}

