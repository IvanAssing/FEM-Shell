#ifndef THINSHELLMESH_H
#define THINSHELLMESH_H

#include "mesh.h"
#include "node.h"
#include "elementsdkt.h"
#include "matrix.h"

class ThinShellMesh : public Mesh
{
    public:
        Node **nodes;
        ElementSDKT **elements;
        Matrix Df,Dm;

        int nNodes;
        int nElements;
        double **results;
        Gnuplot *gnuplot;

        ThinShellMesh(int _nNodes, Node ** _nodes, int _nElements, ElementSDKT **_elements, Matrix _D, Matrix _Dm);

        void plot(vout data);
        virtual void draw(DataGraphic &data);
        virtual void solve(void);

};
#endif // THINSHELLMESH_H
