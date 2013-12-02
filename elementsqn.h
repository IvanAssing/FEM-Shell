#ifndef ELEMENTSQN_H
#define ELEMENTSQN_H

#include "lagrange.h"
#include "node.h"
#include "matrix.h"

class ElementSQN
{
    public:
        int np;
        Node **nodes;
        ElementSQN(int np, Node **nodes);

        void getStiffnessMatrix(Matrix &k, Polynomial2D **Bf, Polynomial2D **Bc, Polynomial2D **Bm, Lagrange &L);

        void draw(void);
};

#endif // ELEMENTSQN_H
