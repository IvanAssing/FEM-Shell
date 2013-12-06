#include "thickplatemesh.h"

#include <fstream>
#include <cstdlib>
#include <QGLWidget>
#include <QDateTime>
#include <cmath>

#define TP_NDOF 3

ThickPlateMesh::ThickPlateMesh(int _nNodes, Node ** _nodes, int _nElements, ElementQN **_elements, int _npx, int _npy, Matrix _D, double _GKt)
    :nNodes(_nNodes), nodes(_nodes), nElements(_nElements), elements(_elements), npx(_npx), npy(_npy), D(_D), GKt(_GKt)
{
    L = new Lagrange(npx-1, npy-1);

    npt = TP_NDOF*npx*npy;

    Polynomial2D Bf[3][npt];
    Polynomial2D Bc[2][npt];

    for(int i=0; i<npx*npy; i++)
    {
        Bf[0][3*i+2] = -1.0*L->D1[i];
        Bf[1][3*i+1] = L->D2[i];
        Bf[2][3*i+1] = L->D1[i];
        Bf[2][3*i+2] = -1.0*L->D2[i];

        Bc[0][3*i] = L->D1[i];
        Bc[0][3*i+2] = L->N[i];
        Bc[1][3*i] = L->D2[i];
        Bc[1][3*i+1] = -1.0*L->N[i];
    }

    BftDBf = new Polynomial2D*[npt];
    for(int i=0; i<npt; i++)
        BftDBf[i] = new Polynomial2D[npt];

    Polynomial2D DB[3][npt];

    for(int i=0; i<3; i++)
        for(int j=0; j<npt; j++)
            DB[i][j] = Bf[0][j]*D(i,0) + Bf[1][j]*D(i,1) + Bf[2][j]*D(i,2);

    for(int i=0; i<npt; i++)
        for(int j=0; j<npt; j++)
            BftDBf[i][j] = Bf[0][i]*DB[0][j] + Bf[1][i]*DB[1][j] + Bf[2][i]*DB[2][j];

    BctBc = new Polynomial2D*[npt];
    for(int i=0; i<npt; i++)
        BctBc[i] = new Polynomial2D[npt];

    for(int i=0; i<npt; i++)
        for(int j=0; j<npt; j++)
            BctBc[i][j] = (Bc[0][i]*Bc[0][j] + Bc[1][i]*Bc[1][j])*GKt;

}



void ThickPlateMesh::solve(void)
{

    int sys_dim = TP_NDOF*nNodes;

    Matrix K(sys_dim , sys_dim );


#pragma omp parallel for
    for(int i=0; i<nElements; i++)
        elements[i]->getStiffnessMatrix(K, BftDBf, BctBc, L);


#pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        if(nodes[i]->lockStatus[2])
            K.setUnit(TP_NDOF*nodes[i]->index + 0);
        if(nodes[i]->lockStatus[3])
            K.setUnit(TP_NDOF*nodes[i]->index + 1);
        if(nodes[i]->lockStatus[4])
            K.setUnit(TP_NDOF*nodes[i]->index + 2);
    }


    Matrix f(sys_dim);


#pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        f(TP_NDOF*nodes[i]->index + 0) = nodes[i]->loadValues[2];
        f(TP_NDOF*nodes[i]->index + 1) = nodes[i]->loadValues[3];
        f(TP_NDOF*nodes[i]->index + 2) = nodes[i]->loadValues[4];
    }

    Matrix x(sys_dim);


    K.solve(f, x);





    Polynomial2D Bf[3][npt];

    Polynomial2D **Bc= new Polynomial2D*[2];
    for(int i=0; i<2; i++)
        Bc[i] = new Polynomial2D[npt];

    for(int i=0; i<npx*npy; i++)
    {
        Bf[0][3*i+2] = -1.0*L->D1[i];
        Bf[1][3*i+1] = L->D2[i];
        Bf[2][3*i+1] = L->D1[i];
        Bf[2][3*i+2] = -1.0*L->D2[i];

        Bc[0][3*i] = GKt * L->D1[i];
        Bc[0][3*i+2] = GKt * L->N[i];
        Bc[1][3*i] = GKt * L->D2[i];
        Bc[1][3*i+1] = -GKt * L->N[i];
    }


    Polynomial2D **DBf= new Polynomial2D*[3];
    for(int i=0; i<npt; i++)
        DBf[i] = new Polynomial2D[npt];


    for(int i=0; i<3; i++)
        for(int j=0; j<npt; j++)
            DBf[i][j] = Bf[0][j]*D(i,0) + Bf[1][j]*D(i,1) + Bf[2][j]*D(i,2);

    Matrix M(sys_dim);
    Matrix Q(2*nNodes);

#pragma omp parallel for
    for(int i=0; i<nElements; i++)
        elements[i]->evalResults(M, Q, x, DBf, Bc);


    results = new double*[2*TP_NDOF+2];
    for(int i=0; i<2*TP_NDOF+2; i++)
        results[i] = new double[nNodes];

#pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        results[0][i] = x(TP_NDOF*i + 0);
        results[1][i] = x(TP_NDOF*i + 1);
        results[2][i] = x(TP_NDOF*i + 2);
        results[3][i] = M(TP_NDOF*i + 0);
        results[4][i] = M(TP_NDOF*i + 1);
        results[5][i] = M(TP_NDOF*i + 2);
        results[6][i] = Q(2*i + 0);
        results[7][i] = Q(2*i + 1);
    }

}


void ThickPlateMesh::draw(DataGraphic &data)
{
    double *x;
    switch (data.var) {
    case W:
        x = results[0];
        break;
    case RX:
        x = results[1];
        break;
    case RY:
        x = results[2];
        break;
    case MX:
        x = results[3];
        break;
    case MY:
        x = results[4];
        break;
    case MXY:
        x = results[5];
        break;
    case QX:
        x = results[6];
        break;
    case QY:
        x = results[7];
        break;
    default:
        return;
        break;
    }

    double T0, T1, T2, T3, T4, Tn, R, G, B;
    getMaxMin(x, nNodes, T4, T0);

    T2 = (T0+T4)/2.;
    T1 = (T0+T2)/2.;
    T3 = (T2+T4)/2.;


    int k[4];

    if(data.def)
        for(int i=0; i<nElements; i++){

            glBegin(GL_QUADS);
            for(int ii = 0; ii<npy-1; ii++)
                for(int jj = 0; jj<npx-1; jj++)
                {
                    k[0] = elements[i]->nodes[ii*npx+jj]->index;
                    k[1] = elements[i]->nodes[ii*npx+jj+1]->index;
                    k[2] = elements[i]->nodes[(ii+1)*npx+jj+1]->index;
                    k[3] = elements[i]->nodes[(ii+1)*npx+jj]->index;

                    for(int p = 0; p<4; p++){
                        Tn = x[k[p]];
                        R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
                        B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
                        G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
                        glColor4d(R,G,B,0.8);

                        glVertex3d(nodes[k[p]]->x, nodes[k[p]]->y, data.factor*results[0][k[p]]);
                    }
                }
            glEnd();

        }
    else

        for(int i=0; i<nElements; i++){

            glBegin(GL_QUADS);
            for(int ii = 0; ii<npy-1; ii++)
                for(int jj = 0; jj<npx-1; jj++)
                {
                    k[0] = elements[i]->nodes[ii*npx+jj]->index;
                    k[1] = elements[i]->nodes[ii*npx+jj+1]->index;
                    k[2] = elements[i]->nodes[(ii+1)*npx+jj+1]->index;
                    k[3] = elements[i]->nodes[(ii+1)*npx+jj]->index;

                    for(int p = 0; p<4; p++){
                        Tn = x[k[p]];
                        R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
                        B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
                        G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
                        glColor4d(R,G,B,0.8);

                        glVertex2d(nodes[k[p]]->x, nodes[k[p]]->y);
                    }
                }
            glEnd();

        }

    if(data.elements)
        for(int i=0; i<nElements; i++)
            elements[i]->draw();

    if(data.load){
        for(int i=0; i<nNodes; i++)
            nodes[i]->draw_lock();

        for(int i=0; i<nNodes; i++)
            nodes[i]->draw_load();
    }

    if(data.nodes){
        for(int i=0; i<nNodes; i++)
            nodes[i]->draw();
    }

}

void ThickPlateMesh::plot(void)
{
    QString zlabel[8];
    zlabel[0] = QString("     w(x,y)");
    zlabel[1] = QString("     Rx(x,y)");
    zlabel[2] = QString("     Ry(x,y)");
    zlabel[3] = QString("     Mx(x,y)");
    zlabel[4] = QString("     My(x,y)");
    zlabel[5] = QString("     Mxy(x,y)");
    zlabel[6] = QString("     Qx(x,y)");
    zlabel[7] = QString("     Qy(x,y)");

    for(int k=0; k<8; k++)
    {
        QDateTime now = QDateTime::currentDateTime();

        QString filename = QString("set output 'graph/FEM-Shell-graphic-")
                + now.toString("yyyyMMddhhmmsszzz") + QString(".png'");

        QString dataname = QString("FEM-Shell-data-")
                + now.toString("yyyyMMddhhmmsszzz") + QString(".tsv");

        std::ofstream file(dataname.toStdString().c_str(),std::ios::out);

        for(int i=0; i<nNodes; i++)
            file<<std::endl<<nodes[i]->x<<"\t"<<nodes[i]->y<<"\t"<<results[k][i];

        file.close();

        Gnuplot g2("points");

        g2.cmd("set terminal pngcairo size 1024,800 enhanced font 'Verdana,10'");

        g2.set_style("points palette pointsize 1 pointtype 7");

        g2.cmd(filename.toStdString());

        g2.set_title("FEM-Shell - Thick Plate Solver");
        g2.set_xlabel("x");
        g2.set_ylabel("y");
        g2.set_zlabel(zlabel[k].toStdString());

        g2.cmd("set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')");

        g2.plotfile_xyz(dataname.toStdString().c_str());
    }

}
