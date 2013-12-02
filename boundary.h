#ifndef BOUNDARY_H
#define BOUNDARY_H

class Boundary
{
    public:

        int index;


        double *loadValues;
        bool *lockStatus;

        Boundary(){}

        Boundary(int index, bool lockStatus[6], double loadValues[6]);
};

#endif // BOUNDARY_H
