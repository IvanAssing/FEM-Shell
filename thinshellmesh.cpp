#include "thinshellmesh.h"

#include <fstream>
#include <cstdlib>
#include <QGLWidget>
#include <QDateTime>

#define FS_NDOF 5


ThinShellMesh::ThinShellMesh(int _nNodes, Node ** _nodes, int _nElements, ElementSDKT **_elements, Matrix _Df, Matrix _Dm)
    :nNodes(_nNodes), nodes(_nodes), nElements(_nElements), elements(_elements), Df(_Df), Dm(_Dm)
{
    gnuplot = new Gnuplot("lines");

}

void ThinShellMesh::solve(void)
{
#pragma omp parallel for num_threads(FEM_NUM_THREADS)
    for(int i=0; i<nElements; i++)
        elements[i]->evaluateTransformationMatrix();

    int sys_dim = FS_NDOF*nNodes;

    Matrix K(sys_dim , sys_dim );

#pragma omp parallel for num_threads(FEM_NUM_THREADS)
    for(int i=0; i<nElements; i++)
        elements[i]->getStiffnessMatrix(K, Df, Dm);

#pragma omp parallel for num_threads(FEM_NUM_THREADS)
    for(int i=0; i<nNodes; i++)
    {
        if(nodes[i]->lockStatus[0])
            K.setUnit(FS_NDOF*nodes[i]->index + 0);
        if(nodes[i]->lockStatus[1])
            K.setUnit(FS_NDOF*nodes[i]->index + 1);
        if(nodes[i]->lockStatus[2])
            K.setUnit(FS_NDOF*nodes[i]->index + 2);
        if(nodes[i]->lockStatus[3])
            K.setUnit(FS_NDOF*nodes[i]->index + 3);
        if(nodes[i]->lockStatus[4])
            K.setUnit(FS_NDOF*nodes[i]->index + 4);
    }

    Matrix f(sys_dim);

#pragma omp parallel for num_threads(FEM_NUM_THREADS)
    for(int i=0; i<nNodes; i++)
    {
        f(FS_NDOF*nodes[i]->index + 0) = nodes[i]->loadValues[0];
        f(FS_NDOF*nodes[i]->index + 1) = nodes[i]->loadValues[1];
        f(FS_NDOF*nodes[i]->index + 2) = nodes[i]->loadValues[2];
        f(FS_NDOF*nodes[i]->index + 3) = nodes[i]->loadValues[3];
        f(FS_NDOF*nodes[i]->index + 4) = nodes[i]->loadValues[4];
    }

    Matrix x(sys_dim);

    K.solve(f, x);

    results = new double*[2*FS_NDOF];
    for(int i=0; i<2*FS_NDOF; i++)
        results[i] = new double[nNodes];

//    Matrix M(sys_dim);

//#pragma omp parallel for num_threads(FEM_NUM_THREADS)
//    for(int i=0; i<nElements; i++)
//        elements[i]->evalResults(M, x, D);


#pragma omp parallel for num_threads(FEM_NUM_THREADS)
    for(int i=0; i<nNodes; i++)
    {
        results[0][i] = x(FS_NDOF*i + 0);
        results[1][i] = x(FS_NDOF*i + 1);
        results[2][i] = x(FS_NDOF*i + 2);
        results[3][i] = x(FS_NDOF*i + 3);
        results[4][i] = x(FS_NDOF*i + 4);
//        results[3][i] = M(FS_NDOF*i + 0);
//        results[4][i] = M(FS_NDOF*i + 1);
//        results[5][i] = M(FS_NDOF*i + 2);
    }

}


void ThinShellMesh::draw(DataGraphic &data)
{

    double *x;
    switch (data.var) {
    case U:
        x = results[0];
        break;
    case V:
        x = results[1];
        break;
    case W:
        x = results[2];
        break;
    case RX:
        x = results[3];
        break;
    case RY:
        x = results[4];
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

    int k[3];

    if(data.def)
        for(int i=0; i<nElements; i++){

            k[0] = elements[i]->n1->index;
            k[1] = elements[i]->n2->index;
            k[2] = elements[i]->n3->index;

            glBegin(GL_TRIANGLES);
            for(int p = 0; p<3; p++){
                Tn = x[k[p]];
                R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
                B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
                G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
                glColor4d(R,G,B,0.8);
                glVertex3d(nodes[k[p]]->x, nodes[k[p]]->y, data.factor*results[0][k[p]]);
            }
            glEnd();

        }
    else
        for(int i=0; i<nElements; i++){

            k[0] = elements[i]->n1->index;
            k[1] = elements[i]->n2->index;
            k[2] = elements[i]->n3->index;

            glBegin(GL_TRIANGLES);
            for(int p = 0; p<3; p++){
                Tn = x[k[p]];
                R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
                B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
                G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
                glColor4d(R,G,B,0.8);
                glVertex2d(nodes[k[p]]->x, nodes[k[p]]->y);
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


void ThinShellMesh::plot(vout data)
{
    int option;
    switch (data) {
    case U:
        option = 0;
        break;
    case V:
        option = 1;
        break;
    case W:
        option = 2;
        break;
    case RX:
        option = 3;
        break;
    case RY:
        option = 4;
        break;
    default:
        return;
        break;
    }


    QString zlabel[5];
    zlabel[0] = QString("u(x,y)");
    zlabel[1] = QString("v(x,y)");
    zlabel[2] = QString("w(x,y)");
    zlabel[3] = QString("Rx(x,y)");
    zlabel[4] = QString("Ry(x,y)");

        QDateTime now = QDateTime::currentDateTime();

        QString dataname = QString("FEM-Shell-data-")
                + now.toString("yyyyMMddhhmmsszzz") + QString(".tsv");

        std::ofstream file(dataname.toStdString().c_str(),std::ios::out);

        for(int i=0; i<nNodes; i++)
            file<<std::endl<<nodes[i]->x<<"\t"<<nodes[i]->y<<"\t"<<results[option][i];

        file.close();


        gnuplot->reset_plot();
        gnuplot->set_style("lines");
        gnuplot->cmd("set dgrid3d 30,30, splines");
        gnuplot->cmd("set hidden3d back offset 1 trianglepattern 3 undefined 1 altdiagonal bentover");
        gnuplot->set_title("FEM-Shell - Thin Flat Shell Solver");
        gnuplot->set_xlabel("x");
        gnuplot->set_ylabel("y");
        gnuplot->set_zlabel(zlabel[option].toStdString());
        gnuplot->set_samples(20);
        gnuplot->set_isosamples(21);
        gnuplot->set_contour();
        gnuplot->unset_legend();
        gnuplot->cmd("set pm3d depthorder");
        gnuplot->cmd("set cntrparam levels auto 20");
        gnuplot->cmd("set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')");
        gnuplot->cmd(QString("splot '%1' u 1:2:3 with pm3d palette").arg(dataname).toStdString());

}
