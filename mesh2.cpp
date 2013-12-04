#include "mesh2.h"

#include "integralgauss2d.h"

#include <fstream>
#include <cstdlib>
#include <QGLWidget>

#include <cmath>

#include <time.h>

Mesh2::Mesh2()
{

//    //clock_t t_start = clock();

//    int np = 2+1;
//    //    int ne = 5;

//    //    int nx = ne*(np-1)+1, ny = ne*(np-1)+1;
//    int ne =10;

//    int nx = ne*(np-1)+1, ny = ne*(np-1)+1;

//    nNodes = nx*ny;

//    nodes = new Node*[nNodes];

//    double dx = 1.0/(nx-1);
//    double dy = 1.0/(ny-1);


//    for(int i=0; i<ny; i++)
//        for(int j=0; j<nx; j++)
//            nodes[j+nx*i] = new Node(j+nx*i, j*dx, i*dy);

//    Node **ptrNodes = new Node*[np*np];


//    nElements = ne*ne;
//    elements = new ElementQN*[nElements];
//    int elementIndex = 0;

//    int ni = 0;
//    for(int ie=0; ie<ne; ie++)
//        for(int je=0; je<ne; je++)
//        {
//            ni = 0;
//            for(int i=0; i<np; i++)
//                for(int j=0; j<np; j++)
//                    ptrNodes[ni++] = nodes[ie*(np-1)*nx + i*nx + je*(np-1) + j];
//            elements[elementIndex++] = new ElementQN(np*np, ptrNodes);
//        }

//    Lagrange L(np-1, np-1);

//    //    std::cout<<L.N[15](1, 1);

//    Polynomial2D Bf[3][3*np*np];
//    Polynomial2D Bc[2][3*np*np];

//    for(int i=0; i<np*np; i++)
//    {
//        Bf[0][3*i+2] = -1.0*L.D1[i];
//        Bf[1][3*i+1] = L.D2[i];
//        Bf[2][3*i+1] = L.D1[i];
//        Bf[2][3*i+2] = -1.0*L.D2[i];

//        Bc[0][3*i] = L.D1[i];
//        Bc[0][3*i+2] = L.N[i];
//        Bc[1][3*i] = L.D2[i];
//        Bc[1][3*i+1] = -1.0*L.N[i];
//    }



//    double vi = 0.3;
//    double E = 200e9;
//    double t = 0.02;
//    double G = 75e9;

//    double GKt = G*5./6.*t;

//    double Ept = E*t*t*t/(12.*(1.0-vi*vi));

//    Matrix D(3,3);

//    D(0, 0) = Ept;
//    D(0, 1) = Ept*vi;
//    D(1, 0) = Ept*vi;
//    D(1, 1) = Ept;
//    D(2, 2) = Ept*(1-vi)/2.0;



//    //    Polynomial2D BftDBf[3*16][3*16];

//    Polynomial2D **BftDBf = new Polynomial2D*[3*np*np];
//    for(int i=0; i<3*np*np; i++)
//        BftDBf[i] = new Polynomial2D[3*np*np];

//    Polynomial2D DB[3][3*np*np];

//    for(int i=0; i<3; i++)
//        for(int j=0; j<3*np*np; j++)
//            DB[i][j] = Bf[0][j]*D(i,0) + Bf[1][j]*D(i,1) + Bf[2][j]*D(i,2);

//    for(int i=0; i<3*np*np; i++)
//        for(int j=0; j<3*np*np; j++)
//            BftDBf[i][j] = Bf[0][i]*DB[0][j] + Bf[1][i]*DB[1][j] + Bf[2][i]*DB[2][j];

//    //    Polynomial2D BctBc[3*16][3*16];
//    Polynomial2D **BctBc = new Polynomial2D*[3*np*np];
//    for(int i=0; i<3*np*np; i++)
//        BctBc[i] = new Polynomial2D[3*np*np];

//    for(int i=0; i<3*np*np; i++)
//        for(int j=0; j<3*np*np; j++)
//            BctBc[i][j] = (Bc[0][i]*Bc[0][j] + Bc[1][i]*Bc[1][j])*GKt;


//    Matrix K(3*nNodes, 3*nNodes);

//    Matrix f(3*nNodes, 1);
//    Matrix x(3*nNodes, 1);

//    #pragma omp parallel for
//        for(int i=0; i<nElements; i++)
//            elements[i]->getStiffnessMatrix(K, BftDBf, BctBc, &L);

//    for(int i=0; i<ny; i++){
//        for(int j=0; j<3*nNodes; j++)
//        {
//            K(3*nodes[i*nx]->index+0, j, 0.0);
//            K(3*nodes[i*nx]->index+1, j, 0.0);
//            K(3*nodes[i*nx]->index+2, j, 0.0);
//        }
//        K(3*nodes[i*nx]->index+0, 3*nodes[i*nx]->index+0, 1.0);
//        K(3*nodes[i*nx]->index+1, 3*nodes[i*nx]->index+1, 1.0);
//        K(3*nodes[i*nx]->index+2, 3*nodes[i*nx]->index+2, 1.0);
//    }

////    for(int i=1; i<=ny; i++)
////        f(3*(i*nx-1), 0, -100.0/ny);

//        f(3*nNodes-3, 0, -1200.);

//        std::cout<<"\n\n *********** SOLVER LINEAR SYSTEM ***********"<<std::flush;

//        std::cout<<"\n\n dim = "<<K.n<<std::flush;

//    K.solve(f, x);

//    std::cout<<"\n\n *********** SOLVER LINEAR SYSTEM ***********"<<std::flush;

//    w = new double[nNodes];
//    for(int i=0; i<nNodes; i++)
//        w[i] = x(3*i, 0);

//    //std::cout<<"\n\n TEMPO = "<<double(clock() - t_start)/CLOCKS_PER_SEC;

//    for(int i=0; i<nNodes; i++)
//        std::cout<<"\n"<<i<<"\t"<<x(3*i+0, 0)<<"\t"<<x(3*i+1, 0)<<"\t"<<x(3*i+2, 0);

//    //plot(x);

}



void Mesh2::plot(Matrix &f)
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

void getMaxMin2(double *vector, int size, double &max, double &min)
{
    max = vector[0];
    min = vector[0];

    for(int i=1; i<size; i++){
        max = vector[i]>max ? vector[i] : max;
        min = vector[i]<min ? vector[i] : min;
    }
}

void Mesh2::draw(void)
{
    double T0, T1, T2, T3, T4, Tn, R, G, B;
    getMaxMin2(w, nNodes, T4, T0);

    T2 = (T0+T4)/2.;
    T1 = (T0+T2)/2.;
    T3 = (T2+T4)/2.;


    int k[4];

    for(int i=0; i<nElements; i++){

        int n = sqrt(elements[i]->np);


        k[0] = elements[i]->nodes[0]->index;
        k[1] = elements[i]->nodes[n-1]->index;
        k[2] = elements[i]->nodes[elements[i]->np-1]->index;
        k[3] = elements[i]->nodes[elements[i]->np-n]->index;

        glBegin(GL_QUADS);
        for(int p = 0; p<4; p++){
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

