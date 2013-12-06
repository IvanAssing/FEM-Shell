#ifndef ELEMENTSSQN_H
#define ELEMENTSSQN_H

#include "lagrange.h"
#include "node.h"
#include "matrix.h"

class ElementSSQN
{
    public:
        int np;
        int npx, npy;
        bool selectiveIntegration;

        Node **nodes;
        ElementSSQN(int npx_, int npy_, Node **nodes_, bool _selectiveIntegracion = false);

        void getStiffnessMatrix(Matrix &k, Polynomial2D **Bf, Polynomial2D **Bc, Polynomial2D **Bm, Lagrange *L);

        void draw(void);
};

#endif // ELEMENTSSQN_H
