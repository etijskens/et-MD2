#include "lj_potential.hpp"

void
compute_forces
  ( py::array_t<double> r
  , py::array_t<double> a
  , py::array_t<double> m
  , VL const & vl
  )
{
    ArrayInfo<double,2> ar(r);
    ArrayInfo<double,2> aa(a);
    ArrayInfo<double,1> am(m);
    std::size_t const n = ar.shape(0);
    double const * ri;
    double const * rj;
    double rij[3];
    for (std::size_t i=0; i<n; ++i)
    {
        ri = &ar[3*i];
        std::size_t const nn = vl.ncontacts(i);
        for (std::size_t ic=0; ic<nn; ++ic)
        {
            std::size_t const j = vl.contact(i,ic);
            rj = &ar[3*j];
            double rij2 = 0;
            for (std::size_t d=0; d<3; ++d) {
                rij[d] = rj[d] - ri[d];
                rij2 += rij[d]*rij[d];
            }
            double ff = lj_force_factor(rij2);
            for (std::size_t d=0; d<3; ++d) {
                aa[3*i+d] += ff*rij[d];
                aa[3*j+d] -= ff*rij[d];
            }
        }
    }
 // convert forces to acceleration
    if (m.shape(0) == 1) {
        for (std::size_t i=0; i<n; ++i) {
            for (std::size_t d=0; d<3; ++d) {
                aa[3*i+d] /= am[0];
            }
        }
    } else if (m.shape(0) == n) {
        for (std::size_t i=0; i<n; ++i) {
            for (std::size_t d=0; d<3; ++d) {
                aa[3*i+d] /= am[i];
            }
        }
    }
}


double
compute_interactions
  ( py::array_t<double> r
  , VL const & vl
  )
{
    ArrayInfo<double,2> ar(r);
    std::size_t const n = ar.shape(0);
    double const * ri;
    double const * rj;
    double epot = 0;
    for (std::size_t i=0; i<n; ++i)
    {
        ri = &ar[3*i];
        std::size_t const nn = vl.ncontacts(i);
        for (std::size_t ic=0; ic<nn; ++ic)
        {
            std::size_t const j = vl.contact(i,ic);
            rj = &ar[3*j];
            double rij2 = 0;
            for (std::size_t d=0; d<3; ++d) {
                rij2 += (rj[d] - ri[d])*(rj[d] - ri[d]);
            }
            epot += lj_potential(rij2);
//            std::cout<<"rij2="<<rij2<<" epot="<<epot<<std::endl;
        }
    }
    return epot;
}

