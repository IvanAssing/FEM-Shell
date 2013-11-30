#ifndef ELEMENTQN_H
#define ELEMENTQN_H

#include "lagrange.h"
#include "node.h"
#include "matrix.h"

class ElementQN
{
    public:
        int np;
        Node **nodes;
        ElementQN(int np, Node **nodes);

        void getStiffnessMatrix(Matrix &k, Polynomial2D **Bf, Polynomial2D **Bc, Lagrange &L);

        void draw(void);
};

#endif // ELEMENTQN_H
