#include "thinplatemesh.h"

#include <fstream>
#include <cstdlib>
#include <QGLWidget>

#define TP_NDOF 3

ThinPlateMesh::ThinPlateMesh(int _nNodes, Node ** _nodes, int _nElements, ElementDKT **_elements, Matrix _D)
    :nNodes(_nNodes), nodes(_nodes), nElements(_nElements), elements(_elements), D(_D)
{

//#pragma omp parallel for
    for(int i=0; i<nElements; i++)
        elements[i]->evaluateTransformationMatrix();

    Matrix K(3*nNodes, 3*nNodes);

    //#pragma omp parallel for
    for(int i=0; i<nElements; i++)
        elements[i]->getStiffnessMatrix(K, D);



    Matrix f(3*nNodes);
    Matrix x(3*nNodes);



#pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        if(nodes[i]->lockStatus[2])
        {
            for(int j=0; j<TP_NDOF*nNodes; j++)
                K(TP_NDOF*nodes[i]->index + 0, j) =  0.0;
            K(TP_NDOF*nodes[i]->index + 0, TP_NDOF*nodes[i]->index + 0) =  1.0;
        }
        if(nodes[i]->lockStatus[3])
        {
            for(int j=0; j<TP_NDOF*nNodes; j++)
                K(TP_NDOF*nodes[i]->index + 1, j) =  0.0;
            K(TP_NDOF*nodes[i]->index + 1, TP_NDOF*nodes[i]->index + 1) =  1.0;
        }
        if(nodes[i]->lockStatus[4])
        {
            for(int j=0; j<TP_NDOF*nNodes; j++)
                K(TP_NDOF*nodes[i]->index + 2, j) =  0.0;
            K(TP_NDOF*nodes[i]->index + 2, TP_NDOF*nodes[i]->index + 2) =  1.0;
        }
    }


#pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        f(TP_NDOF*nodes[i]->index + 0) = nodes[i]->loadValues[2];
        f(TP_NDOF*nodes[i]->index + 1) = nodes[i]->loadValues[3];
        f(TP_NDOF*nodes[i]->index + 2) = nodes[i]->loadValues[4];
    }




    K.solve(f, x);





    w = new double[nNodes];
    for(int i=0; i<nNodes; i++)
        w[i] = x(3*i, 0);

    //    std::cout<<"\n\nw("<<nodes[0]->x<<", "<<nodes[0]->y<<" ) = "<<x(3*nodes[0]->index, 0)<<"\n";

    //    nn = nx*ny/2;
    //    std::cout<<"\n\nw("<<nodes[nn]->x<<", "<<nodes[nn]->y<<" ) = "<<x(3*nodes[nn]->index, 0)<<"\n";
    //std::cout<<"\n\nw("<<nodes[nx*ny-1]->x<<", "<<nodes[nx*ny-1]->y<<" ) = "<<x(3*nodes[nx*ny-1]->index, 0)<<"\n";
    // plot(x);


    //    for(int i=1; i<=ny; i++)
    //        std::cout<<"\n"<<x(3*(i*nx-1), 0);


    //    std::cout<<std::flush;

    //    Node a(0, 0., 0.);
    //    Node c(0, 0., 1.);
    //    Node b(0, 1., 0.);

    //        Node a(0, 0., 0.);
    //        Node b(0, 4., 2.);
    //        Node c(0, 3., 5.);

    //        ElementDKT e(0, &a, &b, &c);

    //        e.evaluateTransformationMatrix();

    //    Matrix ta(3,3), fa(3,1), xa(3,1);

    //    ta(0, 0, 10.);
    //    ta(0, 1, 2.);
    //    ta(1, 0, 3.);
    //    ta(2, 2, 5.);



    //    fa(0, 0, 12.);
    //    fa(1, 0, 3.);
    //    fa(2, 0, 5.);

    //    std::cout<<"\n"<<ta;
    //    std::cout<<"\n"<<fa;
    //    std::cout<<"\n"<<xa;


    //    ta.solve(fa, xa);

    //    std::cout<<"\n\n\n"<<ta;
    //    std::cout<<"\n"<<fa;
    //    std::cout<<"\n"<<xa;





}

void getMaxMin(double *vector, int size, double &max, double &min)
{
    max = vector[0];
    min = vector[0];

    for(int i=1; i<size; i++){
        max = vector[i]>max ? vector[i] : max;
        min = vector[i]<min ? vector[i] : min;
    }
}

void ThinPlateMesh::draw(void)
{
    double T0, T1, T2, T3, T4, Tn, R, G, B;
    getMaxMin(w, nNodes, T4, T0);

    T2 = (T0+T4)/2.;
    T1 = (T0+T2)/2.;
    T3 = (T2+T4)/2.;


    int k[3];

    for(int i=0; i<nElements; i++){

        k[0] = elements[i]->n1->index;
        k[1] = elements[i]->n2->index;
        k[2] = elements[i]->n3->index;

        glBegin(GL_TRIANGLES);
        for(int p = 0; p<3; p++){
            Tn = w[k[p]];
            R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
            B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
            G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
            glColor4d(R,G,B,0.8);

            glVertex2d(nodes[k[p]]->x, nodes[k[p]]->y);
        }
        glEnd();

    }


}


void ThinPlateMesh::plot(Matrix &f)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "plot.png";
    const std::string dat1_filename = "data1.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(int i=0; i<nNodes; i++)
        file1<<nodes[i]->x<<"\t"<<nodes[i]->y<<"\t"<<f(3*i, 0)<<std::endl;
    file1.close();



    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "#set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "#set terminal wxt font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set title \" MESH \"\n"
             "set lmargin 8\n"
             "set xlabel 'x'\n"
             "set ylabel 'y'\n"
             "set zlabel 'w( x, y)'\n"
             "splot '" <<dat1_filename<<"' t\"\" with points lt 2 lc 2 lw 2 pt 3\n"
             "pause -1";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    //const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    //std::system(cmd2.c_str());
}
