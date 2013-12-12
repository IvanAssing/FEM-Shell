#ifndef THICKSHELLMESH_H
#define THICKSHELLMESH_H

#include "mesh.h"
#include "node.h"
#include "elementsqn.h"
#include "lagrange.h"
#include "polynomial2d.h"

class ThickShellMesh : public Mesh
{
    public:
        Lagrange *L;
        Polynomial2D **BftDBf;
        Polynomial2D **BctBc;
        Polynomial2D **BmtDBm;
        int npx, npy;

        Node **nodes;
        ElementSQN **elements;
        Matrix Df, Dm;
        double GKt;
        double **results;
        int nNodes;
        int nElements;
        Gnuplot *gnuplot;

        ThickShellMesh(int _nNodes, Node ** _nodes, int _nElements, ElementSQN **_elements, int _npx, int _npy, Matrix _Df, double _GKt, Matrix _Dm);

        void plot(vout data);
        virtual void draw(DataGraphic &data);
        virtual void solve(void);

    private:
        int npt;
};

#endif // THICKSHELLMESH_H
