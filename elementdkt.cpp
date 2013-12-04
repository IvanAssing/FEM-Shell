#include "elementdkt.h"

#include <QGLWidget>
#include <QFont>
#include <QString>

#define QUADRATIC_TRIANGLE_ { \
    1.0, -3.0, 2.0, -3.0, 4.0, 0.0, 2.0, 0.0, 0.0, \
    0.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
    0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, \
    0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, \
    0.0, 0.0, 0.0, 4.0, -4.0, 0.0, -4.0, 0.0, 0.0, \
    0.0, 4.0, -4.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0}
// 1, x, x², y, yx, yx², y², xy², x²y²


ElementDKT::ElementDKT(int index_, Node *node1, Node *node2, Node *node3)
    :index(index_), n1(node1), n2(node2), n3(node3)
{

}

void ElementDKT::initShapeFunctions(void)
{
    N = new Polynomial2D[6];

    double an[6][9] = QUADRATIC_TRIANGLE_;

    for(int i=0; i<6; i++)
        N[i] = Polynomial2D(2,an[i]);
}

void ElementDKT::evaluateTransformationMatrix(void)
{


    double x[3], y[3], l2[3], a[3], b[3], c[3], d[3], e[3];

    x[0] = n2->x - n3->x;
    x[1] = n3->x - n1->x;
    x[2] = n1->x - n2->x;

    y[0] = n2->y - n3->y;
    y[1] = n3->y - n1->y;
    y[2] = n1->y - n2->y;

    for(int i=0; i<3; i++){
        l2[i] = x[i]*x[i] + y[i]*y[i];
        a[i] = -x[i]/l2[i];
        b[i] = 0.75*x[i]*y[i]/l2[i];
        c[i] = (0.25*x[i]*x[i] - 0.5*y[i]*y[i])/l2[i];
        d[i] = -y[i]/l2[i];
        e[i] = (0.25*y[i]*y[i] - 0.5*x[i]*x[i])/l2[i];
    }

    Polynomial2D Hx[9];
    Polynomial2D Hy[9];

    Hx[0] = 1.5*(a[2]*N[5] - a[1]*N[4]);
    Hx[1] = b[1]*N[4] + b[2]*N[5];
    Hx[2] = N[0] - c[1]*N[4] - c[2]*N[5];

    Hy[0] = 1.5*(d[2]*N[5] - d[1]*N[4]);
    Hy[1] = -1.0*N[0] + e[1]*N[4] + e[2]*N[5];
    Hy[2] = -1.0*Hx[1];

    Hx[3] = 1.5*(a[0]*N[3] - a[2]*N[5]);
    Hx[4] = b[2]*N[5] + b[0]*N[3];
    Hx[5] = N[1] - c[2]*N[5] - c[0]*N[3];

    Hy[3] = 1.5*(d[0]*N[3] - d[2]*N[5]);
    Hy[4] = -1.0*N[1] + e[2]*N[5] + e[0]*N[3];
    Hy[5] = -1.0*Hx[4];

    Hx[6] = 1.5*(a[1]*N[4] - a[0]*N[3]);
    Hx[7] = b[0]*N[3] + b[1]*N[4];
    Hx[8] = N[2] - c[0]*N[3] - c[1]*N[4];

    Hy[6] = 1.5*(d[1]*N[4] - d[0]*N[3]);
    Hy[7] = -1.0*N[2] + e[0]*N[3] + e[1]*N[4];
    Hy[8] = -1.0*Hx[7];

    double by2a = 1.0/(x[1]*y[2] - x[2]*y[1]);

    Polynomial2D dHxd1[9];
    Polynomial2D dHxd2[9];
    Polynomial2D dHyd1[9];
    Polynomial2D dHyd2[9];

    for(int i=0; i<9; i++){
        dHxd1[i] = Hx[i].differential1(1);
        dHxd2[i] = Hx[i].differential2(1);
        dHyd1[i] = Hy[i].differential1(1);
        dHyd2[i] = Hy[i].differential2(1);
    }

    B = new Polynomial2D*[3];
    for(int i=0; i<3; i++)
        B[i] = new Polynomial2D[9];

    for(int i=0; i<9; i++){
        B[0][i] = dHxd1[i]*y[1]*by2a + dHxd2[i]*y[2]*by2a;
        B[1][i] = dHyd1[i]*(-x[1])*by2a + dHyd2[i]*(-x[2])*by2a;
        B[2][i] = dHxd1[i]*(-x[1])*by2a + dHxd2[i]*(-x[2])*by2a + dHyd1[i]*y[1]*by2a + dHyd2[i]*y[2]*by2a;
    }
}


void ElementDKT::getStiffnessMatrix(Matrix &k, Matrix &D)
{
    Polynomial2D BtDB[9][9];
    Polynomial2D DB[3][9];

    for(int i=0; i<3; i++)
        for(int j=0; j<9; j++)
            DB[i][j] = B[0][j]*D(i,0) + B[1][j]*D(i,1) + B[2][j]*D(i,2);

    for(int i=0; i<9; i++)
        for(int j=0; j<9; j++)
            BtDB[i][j] = B[0][i]*DB[0][j] + B[1][i]*DB[1][j] + B[2][i]*DB[2][j];

    // Triangle Gauss integration by 3 points

    int index[3] = {n1->index, n2->index, n3->index};

    // Jabobian/3
    double _2A_by3 = 0.5*((n3->x - n1->x)*(n1->y - n2->y) - (n1->x - n2->x)*(n3->y - n1->y))/3.0; // = 2*A/3

    for(int ii=0; ii<3; ii++)
        for(int ij=0; ij<3; ij++)
            for(int i=0; i<3; i++)
                for(int j=0; j<3; j++)
                    k(3*index[ii] + i, 3*index[ij] + j) += (BtDB[3*ii+i][3*ij+j](0.5, 0.0) +
                            BtDB[3*ii+i][3*ij+j](0.0, 0.5) +
                            BtDB[3*ii+i][3*ij+j](0.5, 0.5))*_2A_by3;
}

void ElementDKT::draw(void)
{

    glColor4d(0.0, 1.0, 0.0, 0.8);
    glLineWidth(2.5f);

    glBegin(GL_LINE_LOOP);
    {
        glVertex3d(n1->x, n1->y, 0.0);
        glVertex3d(n2->x, n2->y, 0.0);
        glVertex3d(n3->x, n3->y, 0.0);
    }
    glEnd();
}
