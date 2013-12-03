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
        //double operator()(int i, int j);
        void operator = (Matrix &M);

        double& operator()(int i, int j)
        {
          return data[j*n + i];
        }

        const double& operator()(int i, int j) const
        {
          return data[j*n + i];
        }


        void solve(Matrix &b, Matrix &x);

        friend std::ostream& operator<<(std::ostream& out, Matrix& obj);

        static void evalInverseMatrix(int n_, double **A);

private:
        int dim;
        double *data;
};

#endif // MATRIX_H
