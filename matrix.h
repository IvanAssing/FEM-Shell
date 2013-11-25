#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

class Matrix
{
    public:
        int m,n;

        Matrix();
        Matrix(int m, int n);
        void operator()(int i, int j, double value);
        void add(int i, int j, double value);
        double operator()(int i, int j);


        void solve(Matrix &b, Matrix &x);

        friend std::ostream& operator<<(std::ostream& out, Matrix& obj);

        static void evalInverseMatrix(int n_, double **A);

private:
        int dim;
        double *data;
};

#endif // MATRIX_H
