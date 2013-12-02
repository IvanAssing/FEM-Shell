#include "boundary.h"

Boundary::Boundary(int index_, bool lockStatus_[6], double loadValues_[6])
:index(index_)
{
    loadValues = new double[6];
    lockStatus = new bool[6];

    for(int i=0; i<6; i++)
    {
        loadValues[i] = loadValues_[i];
        lockStatus[i] = lockStatus_[i];
    }
}
