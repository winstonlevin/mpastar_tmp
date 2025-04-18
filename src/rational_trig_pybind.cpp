#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "rational_trig.h"

namespace py = pybind11;

void init_rational_trig(py::module_ m) {
    py::class_<RationalTrig>(m, "RationalTrig")
        .def(py::init())
        .def(py::init<unsigned short int>())
        .def(py::init<unsigned short int, std::vector<unsigned short int>, std::vector<unsigned short int>, std::vector<unsigned short int>>())
        .def(py::init<unsigned short int, std::vector<unsigned short int>, std::vector<unsigned short int>, std::vector<unsigned short int>, std::vector<double>>())
        
        .def_readwrite("den", &RationalTrig::den)
        .def_readwrite("ns", &RationalTrig::ns)
        .def_readwrite("nc", &RationalTrig::nc)
        .def_readwrite("fac", &RationalTrig::fac)
        .def_readwrite("tha", &RationalTrig::tha)
        .def("__repr__", &RationalTrig::to_string)
        ;
}  // init_rational_trig
