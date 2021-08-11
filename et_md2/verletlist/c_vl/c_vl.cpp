/*
 *  C++ source file for module et_md2.verletlist.c_vl
 */


#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include "vl_lib/vl.cpp"

PYBIND11_MODULE(c_vl, m)
{// doc-string
    m.doc() = "Verlet List, C++ implementation";

 // exposed classes
 // Verlet list
    py::class_<VL>(m, "VL")
        .def(py::init<std::size_t, double>())
        .def("reset"     , &VL::reset)
        .def("add"       , &VL::add)
        .def("linearise" , &VL::linearise)
        .def("natoms"    , &VL::natoms)
        .def("has"       , &VL::has)
        .def("print"     , &VL::print)
        .def("contact"   , &VL::contact)
        .def("ncontacts" , &VL::ncontacts)
        .def("cutoff"   , &VL::cutoff)
    ;
 // Hilbert curve functions
//    m.def("xyzw2h_float64", xyzw2h_float64 ); // 3D positions to hilbert index of the corresponding cell with width w.
//    m.def("xyzw2h_float32", xyzw2h_float32 ); // 3D positions to hilbert index of the corresponding cell with width w.
//
//    m.def("xyzw2ijkh_float32", xyzw2ijkh_float32 ); // 3D positions to cell indices and hilbert index of the corresponding cell with width w.
//    m.def("xyzw2ijkh_float64", xyzw2ijkh_float64 ); // 3D positions to cell indices and hilbert index of the corresponding cell with width w.
//
//    m.def(  "rw2h_float64",   rw2h_float64 ); // 3D positions to hilbert index of the corresponding cell with width w.
//    m.def(  "rw2h_float32",   rw2h_float32 ); // 3D positions to hilbert index of the corresponding cell with width w.
//
//    m.def("rw2ch_float64", rw2ch_float64 ); // 3D positions to cell indices and hilbert index of the corresponding cell with width w.
//    m.def("rw2ch_float32", rw2ch_float32 ); // 3D positions to cell indices and hilbert index of the corresponding cell with width w.

 // spatial sorting
//    m.def("sort"  , sort   ); // sort the hilbert indices and produce reordering array.
//    m.def("reorder_float32"     , reorder_float32);
//    m.def("reorder_float64"     , reorder_float64);
//    m.def("reorder_int32"       , reorder_int32  );
//    m.def("reorder_uint32"      , reorder_uint32 );
//    m.def("reorder_longlongint" , reorder_longlongint );
//    m.def("reorder2_float64"    , reorder2_float64);
//
//    m.def("ijk2h_1",ijk2h_1);
//    m.def("h2ijk_1",h2ijk_1);
//    m.def("h2ijk",h2ijk);
//
//    m.def("info",hilbert::info);
//    m.def("cell_index_limit",hilbert::cell_index_limit);
//    m.def("hilbert_index_limit",hilbert::hilbert_index_limit);
//    m.def("validate_cell_index",hilbert::validate_cell_index);
//    m.def("validate_hilbert_index",hilbert::validate_hilbert_index);
//    m.def("is_validating",hilbert::is_validating);
//    m.def("is_valid_ijk",hilbert::is_valid_ijk);
//    m.def("is_valid_i",hilbert::is_valid_i);
//    m.def("is_valid_h",hilbert::is_valid_h);
//
//    m.def("build_hl", build_hl);
//    m.def("build_vl", build_vl);
}

