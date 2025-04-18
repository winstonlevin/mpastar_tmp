#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


namespace py = pybind11;

void init_dubins(py::module_ m);
void init_astar(py::module_ m);

PYBIND11_MODULE(_core, m) {
    m.doc() = "mpastar";

    py::module_ dubins = m.def_submodule("dubins", "Dubins submodule of mpastar");
    init_dubins(dubins);

    py::module_ graphs = m.def_submodule("graphs", "Graph search submodule of mpastar");
    init_astar(graphs);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
