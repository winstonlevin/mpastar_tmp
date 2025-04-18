#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "astar.h"

namespace py = pybind11;

void init_astar(py::module_ m) {
    py::class_<IntegerPosition>(m, "IntegerPosition")
        .def(py::init())
        .def(py::init<int, int>())
        .def_readwrite("x", &IntegerPosition::x)
        .def_readwrite("y", &IntegerPosition::y)
        .def("__repr__", &IntegerPosition::to_string)
        ;

    py::class_<Node>(m, "Node")
        .def(py::init())
        .def(py::init<int, int>())
        .def(py::init<int, int, int, int>())
        .def(py::init<int, int, int, int, int>())
        .def(py::init<int, int, int, int, int, bool>())
        .def(py::init<int, int, int, int, int, bool, bool>())
        .def_readwrite("i", &Node::i)
        .def_readwrite("j", &Node::j)
        .def_readwrite("g", &Node::g)
        .def_readwrite("h", &Node::h)
        .def_readwrite("f", &Node::f)
        .def_readwrite("ip", &Node::ip)
        .def_readwrite("jp", &Node::jp)
        .def_readwrite("open", &Node::open)
        .def_readwrite("terminal", &Node::terminal)
        .def_readwrite("obstructed", &Node::obstructed)
        .def("__repr__", &Node::to_string)
        .def("update_cost_to_come", &Node::update_cost_to_come)
        ;

    py::class_<ManhattanAstar>(m, "ManhattanAstar")
        .def(py::init())
        .def(py::init<IntegerPosition, IntegerPosition>())
        .def(py::init<IntegerPosition, IntegerPosition, int, int>())
        .def(py::init<IntegerPosition, IntegerPosition, int, int, int>())
        .def(py::init<IntegerPosition, IntegerPosition, std::vector<IntegerPosition>>())
        .def(py::init<IntegerPosition, IntegerPosition, std::vector<IntegerPosition>, int, int>())
        .def(py::init<IntegerPosition, IntegerPosition, std::vector<IntegerPosition>, int, int, int>())
        .def_readwrite("p0", &ManhattanAstar::p0)
        .def_readwrite("pf", &ManhattanAstar::pf)
        .def_readwrite("height", &ManhattanAstar::height)
        .def_readwrite("width", &ManhattanAstar::width)
        .def_readwrite("max_iter", &ManhattanAstar::max_iter)
        .def_readwrite("status", &ManhattanAstar::status)
        .def_readwrite("nodes", &ManhattanAstar::nodes)
        .def_readwrite("obstacles", &ManhattanAstar::obstacles)
        .def("__repr__", &ManhattanAstar::to_string)
        .def("solve_graph", &ManhattanAstar::solve_graph, "Apply A* algorithm to fill in cost-to-come until parent to goal node is discovered")
        .def("construct_path", &ManhattanAstar::construct_path, "Construct path from initial node to node i, j")
        ;
    
    py::class_<AstarNode::AstarEdge>(m, "AstarEdge")
        .def(py::init())
        .def(py::init<AstarNode*, float>())
        .def_readwrite("child", &AstarNode::AstarEdge::child)
        .def_readwrite("edge_cost", &AstarNode::AstarEdge::edge_cost)
        .def("__repr__", &AstarNode::AstarEdge::to_string)
        ;

    py::class_<AstarNode>(m, "AstarNode")
        .def(py::init())
        .def(py::init<unsigned short int, float, bool, bool>())
        .def(py::init<unsigned short int, std::vector<AstarNode::AstarEdge>, float, bool, bool>())
        .def_readwrite("id", &AstarNode::id)
        .def_readwrite("parent", &AstarNode::parent)
        .def_readwrite("edges", &AstarNode::edges)
        .def_readwrite("cost_to_come", &AstarNode::cost_to_come)
        .def_readwrite("heuristic", &AstarNode::heuristic)
        .def_readwrite("cost", &AstarNode::cost)
        .def_readwrite("open", &AstarNode::open)
        .def_readwrite("terminal", &AstarNode::terminal)
        .def_readwrite("req_init", &AstarNode::req_init)
        .def("__repr__", &AstarNode::to_string)
        .def("check_parent", &AstarNode::check_parent)
        .def(py::self > py::self)
        .def(py::self < py::self)
        .def(py::self == py::self)
        ;

    m.def("construct_path", &construct_path, "Generate path backward from Node");
}  // init_astar
