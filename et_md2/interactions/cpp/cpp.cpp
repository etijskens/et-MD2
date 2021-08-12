/*
 *  C++ source file for module et_md2.interactions.cpp
 */


// See http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html for examples on how to use pybind11.
// The example below is modified after http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html#More-on-working-with-numpy-arrays
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

#include "vl.hpp"
#include "ArrayInfo.hpp"


template<typename FloatType, typename Potential>
void
compute_forces
  ( py::array_t<FloatType> r
  , py::array_t<FloatType> a
  , py::array_t<FloatType> m
  , VL const & vl
  )
{
    ArrayInfo<FloatType,2> ar(r);
    ArrayInfo<FloatType,2> aa(a);
    ArrayInfo<FloatType,1> am(m);
    std::size_t const n = ar.shape(0);
    FloatType const * ri;
    FloatType const * rj;
    FloatType rij[3];
    for (std::size_t i=0; i<n; ++i)
    {
        ri = &ar[3*i];
        std::size_t const nn = vl.ncontacts(i);
        for (std::size_t ic=0; ic<nn; ++ic)
        {
            std::size_t const j = vl.contact(i,ic);
            rj = &ar[3*j];
            FloatType rij2 = 0;
            for (std::size_t d=0; d<3; ++d) {
                rij[d] = rj[d] - ri[d];
                rij2 += rij[d]*rij[d];
            }
            FloatType ff = force_factor(rij2);
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

template <typename FloatType>
FloatType
compute_interactions
  ( py::array_t<FloatType> r
  , VL const & vl
  )
{
    ArrayInfo<FloatType,2> ar(r);
    std::size_t const n = ar.shape(0);
    FloatType const * ri;
    FloatType const * rj;
    FloatType epot = 0;
    for (std::size_t i=0; i<n; ++i)
    {
        ri = &ar[3*i];
        std::size_t const nn = vl.ncontacts(i);
        for (std::size_t ic=0; ic<nn; ++ic)
        {
            std::size_t const j = vl.contact(i,ic);
            rj = &ar[3*j];
            FloatType rij2 = 0;
            for (std::size_t d=0; d<3; ++d) {
                rij2 += (rj[d] - ri[d])*(rj[d] - ri[d]);
            }
            epot += potential(rij2);
//            std::cout<<"rij2="<<rij2<<" epot="<<epot<<std::endl;
        }
    }
    return epot;
}


PYBIND11_MODULE(cpp, m)
{// optional module doc-string
    m.doc() = "C++ implementation of et_md2.interactions"; // optional module docstring
 // list the functions you want to expose:
    m.def("compute_forces_sp"      , &compute_forces<float>);
    m.def("compute_interactions_sp", &compute_interactions<float>);
    m.def("compute_forces_dp"      , &compute_forces<double>);
    m.def("compute_interactions_dp", &compute_interactions<double>);
}
