#ifndef NODE_H
#define NODE_H

#include "boundary.h"

class Node
{
    public:
        int index;
        double x;
        double y;
        double z;

        double *loadValues;
        bool *lockStatus;


        Node();
        Node(int index, double x, double y, double z = 0.0);

        void setup(Boundary &b);
        void draw(void);
        void draw_lock(void);
        void draw_load(void);

        ~Node();
};

#endif // NODE_H
