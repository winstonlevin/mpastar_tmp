#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "dubins.h"

namespace py = pybind11;

void init_dubins(py::module_ m) {
    py::class_<Position>(m, "Position")
        .def(py::init())
        .def(py::init<double, double>())
        .def_readwrite("x", &Position::x)
        .def_readwrite("y", &Position::y)
        .def(py::self + py::self)
        .def(py::self + double())
        .def(double() + py::self)
        .def(py::self += py::self)
        .def(py::self += double())
        .def(py::self - py::self)
        .def(py::self - double())
        .def(double() - py::self)
        .def(py::self -= py::self)
        .def(py::self -= double())
        .def(-py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self *= double())
        .def(py::self / double())
        .def(double() / py::self)
        .def(py::self /= double())
        .def("__repr__", &Position::to_string)
        .def("norm2", &Position::norm2, "Squared cartesian norm of position")
        .def("norm", &Position::norm, "Cartesian norm of position")
        .def("angle", &Position::angle, "Angle of y/x")
        .def("dot", &Position::dot, "Dot product of 2 positions")
        ;

    py::class_<Angle>(m, "Angle")
        .def(py::init())
        .def(py::init<double>())
        .def(py::init<double, double, double>())
        .def_readwrite("a", &Angle::a)
        .def_readwrite("sina", &Angle::sina)
        .def_readwrite("cosa", &Angle::cosa)
        .def("__repr__", &Angle::to_string)
        .def("tana", &Angle::tana)
        .def("cota", &Angle::cota)
        .def("seca", &Angle::seca)
        .def("csca", &Angle::csca)
        .def("negate", &Angle::negate)
        .def("negate_in_place", &Angle::negate_in_place)
        .def("rotate", &Angle::rotate)
        .def("rotate_in_place", &Angle::rotate_in_place)
        .def("rotate_90", &Angle::rotate_90)
        .def("rotate_180", &Angle::rotate_180)
        .def("rotate_270", &Angle::rotate_270)
        .def("rotate_90_in_place", &Angle::rotate_90_in_place)
        .def("rotate_180_in_place", &Angle::rotate_180_in_place)
        .def("rotate_270_in_place", &Angle::rotate_270_in_place)
        .def(py::self + py::self)
        .def(double() + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(double() - py::self)
        .def(py::self -= py::self)
        .def(-py::self)
        ;

    m.def("position_to_angle", &position_to_angle, "Convert postion holding cos(psi), sin(psi) to angle psi");
    m.def("angle_to_position", &angle_to_position, "Extract position holding cos(psi), sin(psi) from angle psi");

    py::class_<State>(m, "State")
        .def(py::init())
        .def(py::init<Position, Angle>())
        .def(py::init<Position, double>())
        .def(py::init<Position, double, double, double>())
        .def(py::init<double, double, Angle>())
        .def(py::init<double, double, double>())
        .def(py::init<double, double, double, double, double>())
        .def_readwrite("p", &State::p)
        .def_readwrite("psi", &State::psi)
        .def("__repr__", &State::to_string)
        ;

    py::class_<Turn>(m, "Turn")
        .def(py::init())
        .def(py::init<State, double>())
        .def(py::init<State, double, unsigned short int>())
        .def_readwrite("s", &Turn::s)
        .def_readwrite("r", &Turn::r)
        .def_readwrite("nf_id", &Turn::nf_id)
        .def("__repr__", &Turn::to_string)
        .def("pc", &Turn::pc)
        ;

    py::class_<DubinsPath>(m, "DubinsPath")
        .def(py::init())
        .def(py::init<Turn, Turn, bool, std::vector<double>, State, State>())
        .def(py::init<Turn, Turn, bool, std::vector<double>, State, State, double>())
        .def_readwrite("t0", &DubinsPath::t0)
        .def_readwrite("tf", &DubinsPath::tf)
        .def_readwrite("x1", &DubinsPath::x1)
        .def_readwrite("x2", &DubinsPath::x2)
        .def_readwrite("csc", &DubinsPath::csc)
        .def_readwrite("tha", &DubinsPath::tha)
        .def_readwrite("cost", &DubinsPath::cost)
        .def("__repr__", &DubinsPath::to_string)
        .def(py::self > py::self)
        .def(py::self < py::self)
        .def(py::self == py::self)
        ;

    m.def("pc2pw", &pc2pw, "Convert center of turn state to tangency position");
    m.def("pw2pc", &pw2pc, "Convert tangency state to center of turn");
    m.def("law_of_cos", &law_of_cos, "Return the cosine of angle for triangle with sides equal to ``side_opp'', ``side_1'', and ``side_2''. The angle is opposite to ``side_opp''");
    m.def("wrap_dpsi", py::overload_cast<double, double, double>(&wrap_dpsi), "Wrap heading based on turn direction");
    m.def("wrap_dpsi", py::overload_cast<double, double>(&wrap_dpsi), "Wrap heading based on turn direction");
    m.def("wrap_to_pi", &wrap_to_pi, "Wrap angle to be in range [-pi, pi]");
    m.def("connect_configuration", py::overload_cast<State, double, State, double>(&connect_configuration));
    m.def("connect_configuration", py::overload_cast<State, State>(&connect_configuration));
    m.def("connect_csc", py::overload_cast<Turn, Turn, double>(&connect_csc));
    m.def("connect_csc", py::overload_cast<Turn, Turn>(&connect_csc));
    m.def("connect_ccc", py::overload_cast<Turn, Turn, double>(&connect_ccc));
    m.def("connect_ccc", py::overload_cast<Turn, Turn>(&connect_ccc));

    py::class_<LineSegment>(m, "LineSegment")
        .def(py::init())
        .def(py::init<Position, Position>())
        .def_readwrite("p0", &LineSegment::p0)
        .def_readwrite("v", &LineSegment::v)
        .def("__repr__", &LineSegment::to_string)
        ;

    py::class_<Circle>(m, "Circle")
        .def(py::init())
        .def(py::init<Position>())
        .def(py::init<Position, double>())
        .def_readwrite("pc", &Circle::pc)
        .def_readwrite("rc", &Circle::rc)
        .def_readwrite("rc2", &Circle::rc2)
        .def("__repr__", &Circle::to_string)
        ;

    py::class_<Polygon>(m, "Polygon")
        .def(py::init())
        .def(py::init<std::vector<LineSegment>>())
        .def_readwrite("edges", &Polygon::edges)
        .def_readwrite("pc", &Polygon::pc)
        .def_readwrite("rc", &Polygon::rc)
        .def_readwrite("rc2", &Polygon::rc2)
        .def("__repr__", &Polygon::to_string)
        .def("calculate_neighborhood", &Polygon::calculate_neighborhood)
        .def("position_inside", &Polygon::position_inside)
        .def("line_intersect", py::overload_cast<const LineSegment, const double>(&Polygon::line_intersect, py::const_))
        .def("line_intersect", py::overload_cast<const LineSegment>(&Polygon::line_intersect, py::const_))
        .def("turn_intersect", py::overload_cast<const Position, const double, const Position, const Position, const double>(&Polygon::turn_intersect, py::const_))
        .def("turn_intersect", py::overload_cast<const Position, const Position, const Position, const double>(&Polygon::turn_intersect, py::const_))
        .def("turn_intersect", py::overload_cast<const Position, const double, const Position, const Position>(&Polygon::turn_intersect, py::const_))
        .def("turn_intersect", py::overload_cast<const Position, const Position, const Position>(&Polygon::turn_intersect, py::const_))
        .def("convex_circumscription", &Polygon::convex_circumscription)
        ;

    py::class_<ConvexHull>(m, "ConvexHull")
        .def(py::init())
        .def(py::init<std::vector<LineSegment>>())
        .def_readwrite("edges", &ConvexHull::edges)
        .def("__repr__", &ConvexHull::to_string)
        .def("point_intersect", py::overload_cast<const Position, const double>(&ConvexHull::point_intersect, py::const_), "Check if point lies inside convex hull")
        .def("point_intersect", py::overload_cast<const Position>(&ConvexHull::point_intersect, py::const_), "Check if point lies inside convex hull")
        .def("line_intersect", py::overload_cast<const Position, const Position, const double>(&ConvexHull::line_intersect, py::const_), "Check if point lies inside convex hull")
        .def("line_intersect", py::overload_cast<const Position, const Position>(&ConvexHull::line_intersect, py::const_), "Check if point lies inside convex hull")
        ;

    m.def("plane_cut", py::overload_cast<const Position, const LineSegment>(&plane_cut), "Divide point into half-plane cut by edge, negative values indicating point is inside half plane");
    m.def("plane_cut", py::overload_cast<const Position, const Position>(&plane_cut), "Divide point into half-plane cut by edge, negative values indicating point is inside half plane");
    m.def("plane_cut", py::overload_cast<const Position, const Position, const Position>(&plane_cut), "Divide point into half-plane cut by edge, negative values indicating point is inside half plane");
    m.def("trig_in_range_wrapped", py::overload_cast<const Position, const Position, const Position, const double>(&trig_in_range_wrapped), "Determine if wrapped angle is in range using (cos,sin) pairs");
    m.def("trig_in_range_wrapped", py::overload_cast<const Position, const Position, const Position>(&trig_in_range_wrapped), "Determine if wrapped angle is in range using (cos,sin) pairs");
    m.def("trig_in_range_wrapped", py::overload_cast<const Angle, const Position, const Position, const double>(&trig_in_range_wrapped), "Determine if wrapped angle is in range using (cos,sin) pairs");
    m.def("trig_in_range_wrapped", py::overload_cast<const Angle, const Position, const Position>(&trig_in_range_wrapped), "Determine if wrapped angle is in range using (cos,sin) pairs");
}  // init_dubins
