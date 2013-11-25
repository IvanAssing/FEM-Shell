#include "polynomial1d.h"

#include <iomanip>
#include <omp.h>

Polynomial1D::Polynomial1D()
    :n(1)
{
    an = new double[1];
    an[0] = 0.0;
}


Polynomial1D::Polynomial1D(const Polynomial1D& b)
    :n(b.n)
{
    an = new double[b.n];

    //#pragma omp parallel for
    for(int i=0; i<n; i++)
        an[i] = b.an[i];
}


Polynomial1D::Polynomial1D(int degree, double *coefficients)
    :n(degree+1)
{
    an = new double[n];

//#pragma omp parallel for
    for(int i=0; i<n; i++)
        an[i] = coefficients[i];

}


int Polynomial1D::getDegree(void) const
{
    return n-1;
}


double Polynomial1D::integral(const double xinf, const double xsup) const
{
    double xa = xinf;
    double xb = xsup;
    double dn = 1.0;
    double sum = an[0]*(xb-xa);

    for(int i=1; i<n; i++){
        xa *= xinf;
        xb *= xsup;
        dn += 1.0;
        sum += an[i]*(xb-xa)/dn;
    }

    return sum;
}


Polynomial1D& Polynomial1D::differential(int nDiffOrder) const
{
    if(nDiffOrder<1) return *new Polynomial1D(*this);

    Polynomial1D *out = new Polynomial1D;

    if(nDiffOrder>=n) return *out;

    out->resize(this->n - nDiffOrder);

    double fat = 1.0;
    for(int i=0; i<nDiffOrder; i++)
        fat *= static_cast<double>(n - i - 1);

    for(int i=n-1; i-nDiffOrder>=0; i--){
        out->an[i-nDiffOrder] = an[i]*fat;
        fat *= static_cast<double>(i-nDiffOrder)/static_cast<double>(i);
    }

    return *out;

}



void Polynomial1D::resize(int newSize)
{
    this->n = newSize;
    delete [] an;

    an = new double[n];

//#pragma omp parallel for
    for(int i=0; i<n; i++)
        an[i] = 0.0;

}


Polynomial1D::~Polynomial1D()
{
    delete [] an;
}



double Polynomial1D::operator()(double x) const
{
    double sum = an[0];
    double pxn = 1.0;

    for(int i=1; i<n; i++)
    {
        pxn *= x;
        sum += an[i]*pxn;
    }

    return sum;
}


Polynomial1D& Polynomial1D::operator=(const Polynomial1D& px)
{
    this->resize(px.n);

//#pragma omp parallel for
    for(int i=0; i<px.n; i++)
        an[i] = px.an[i];

    return *this;
}


Polynomial1D& Polynomial1D::operator *(const Polynomial1D& b) const
{
    Polynomial1D *prod = new Polynomial1D;

    if(n<2)
    {
        prod->resize(b.n);
//#pragma omp parallel for
        for(int i=0; i<b.n; i++)
            prod->an[i] = b.an[i]*an[0];
        return *prod;
    }

    if(b.n<2)
    {
        prod->resize(n);

        //#pragma omp parallel for
        for(int i=0; i<n; i++)
            prod->an[i] = an[i]*b.an[0];
        return *prod;
    }

    prod->resize(n + b.n - 1);

//#pragma omp parallel for
    for(int i=0; i<n; i++)
        for(int j=0; j<b.n; j++)
            prod->an[i+j] += an[i]*b.an[j];

    return *prod;

}


Polynomial1D& Polynomial1D::operator+(const Polynomial1D& px) const
{
    Polynomial1D *sum = new Polynomial1D;

    if(n>=px.n)
    {
        sum->resize(n);

//#pragma omp parallel for
        for(int i=0; i<n; i++)
            sum->an[i] = an[i];

//#pragma omp parallel for
        for(int i=0; i<px.n; i++)
            sum->an[i] += px.an[i];

        return *sum;
    }
    else
    {
        sum->resize(px.n);

//#pragma omp parallel for
        for(int i=0; i<px.n; i++)
            sum->an[i] = px.an[i];

//#pragma omp parallel for
        for(int i=0; i<n; i++)
            sum->an[i] += an[i];

        return *sum;
    }
}


Polynomial1D& Polynomial1D::operator-(const Polynomial1D& px) const
{
    Polynomial1D *sum = new Polynomial1D;

    if(n>=px.n)
    {
        sum->resize(n);

//#pragma omp parallel for
        for(int i=0; i<n; i++)
            sum->an[i] = an[i];

//#pragma omp parallel for
        for(int i=0; i<px.n; i++)
            sum->an[i] -= px.an[i];

        return *sum;
    }
    else
    {
        sum->resize(px.n);

        //#pragma omp parallel for
        for(int i=0; i<px.n; i++)
            sum->an[i] = -px.an[i];

        //#pragma omp parallel for
        for(int i=0; i<n; i++)
            sum->an[i] += an[i];

        return *sum;
    }
}


Polynomial1D& Polynomial1D::operator*(const double alpha) const
{
    Polynomial1D *prod = new Polynomial1D(*this);

        //#pragma omp parallel for
    for(int i=0; i<prod->n; i++)
        prod->an[i] *= alpha;

    return *prod;
}


namespace MC{

    std::ostream&  operator<<(std::ostream& out, Polynomial1D& px)
    {
        out<<std::setprecision(4)<<std::setw(10)<<std::uppercase;
        out<<std::scientific<<px.an[0];

        for(int i=1; i<px.n; i++)
            out<<" + "<<px.an[i]<<" X^"<<i;

        return out;
    }

}
