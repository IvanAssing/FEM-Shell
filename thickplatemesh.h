#ifndef THICKPLATEMESH_H
#define THICKPLATEMESH_H

#include "mesh.h"
#include "node.h"
#include "elementqn.h"
#include "lagrange.h"
#include "polynomial2d.h"

class ThickPlateMesh : public Mesh
{
    public:
        Lagrange *L;
        Polynomial2D **BftDBf;
        Polynomial2D **BctBc;
        int npx, npy;
        Node **nodes;
        ElementQN **elements;
        Matrix D;
        double GKt;
        double **results;
        Gnuplot *gnuplot;
        int nNodes;
        int nElements;

        ThickPlateMesh(int _nNodes, Node ** _nodes, int _nElements, ElementQN **_elements, int _npx, int _npy, Matrix _D, double _GKt);

        void plot(vout data);
        virtual void draw(DataGraphic &data);
        virtual void solve(void);

    private:
        int npt;
};

#endif // THICKPLATEMESH_H
