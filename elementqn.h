#ifndef ELEMENTQN_H
#define ELEMENTQN_H

#include "lagrange.h"
#include "node.h"
#include "matrix.h"

class ElementQN
{
    public:
        Node **nodes;
        ElementQN(Node **nodes);

        void getStiffnessMatrix(Matrix &k, Polynomial2D **Bf, Polynomial2D **Bc);

        void draw(void);
};

#endif // ELEMENTQN_H
