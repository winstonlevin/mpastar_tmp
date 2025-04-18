#include "dubins.h"

using namespace std::numbers;

const double DUB_TOL = std::numeric_limits<double>::epsilon();
const double DUB_INF = std::numeric_limits<double>::infinity();
bool dub_eq(double a, double b) {
    return fabs(a - b) < DUB_TOL;
}

// C++ fmod does NOT follow same convention as Python's mod operator: rounds toward 0 instead of to same sign as b,
// so ``ufmod'' is used instead to wrap the angle
inline double ufmod(double a, double b) {return fmod(fmod(a, b) + b, b);}

// Configuration structures
Position Position::operator+(const Position& other) const {
    return Position(x + other.x, y + other.y);
}
Position Position::operator+(const double& other) const {
    return Position(x + other, y + other);
}
Position operator+(double other, const Position &p) {
    return Position(other + p.x, other + p.y);
}
Position& Position::operator+=(const Position& other) {
    x += other.x;
    y += other.y;
    return *this;
}
Position& Position::operator+=(const double& other) {
    x += other;
    y += other;
    return *this;
}
Position Position::operator-(const Position& other) const {
    return Position(x - other.x, y - other.y);
}
Position Position::operator-(const double& other) const {
    return Position(x - other, y - other);
}
Position Position::operator-() const {
    return Position(-x, -y);
}
Position operator-(double other, const Position &p) {
    return Position(other - p.x, other - p.y);
}
Position& Position::operator-=(const Position& other) {
    x -= other.x;
    y -= other.y;
    return *this;
}
Position& Position::operator-=(const double& other) {
    x -= other;
    y -= other;
    return *this;
}
Position Position::operator*(const double& other) const {
    return Position(x * other, y * other);
}
Position operator*(double other, const Position &p) {
    return Position(other * p.x, other * p.y);
}
Position& Position::operator*=(const double& other) {
    x *= other;
    y *= other;
    return *this;
}
Position Position::operator/(const double& other) const {
    return Position(x / other, y / other);
}
Position operator/(double other, const Position &p) {
    return Position(other / p.x, other / p.y);
}
Position& Position::operator/=(const double& other) {
    x /= other;
    y /= other;
    return *this;
}

std::string Position::to_string() const {
    return "Position[" + std::to_string(x) + ", " + std::to_string(y) + "]";
}

double Position::norm2() const {
    return x*x + y*y;
}
double Position::norm() const {
    return sqrt(norm2());
}
double Position::angle() const {
    return atan2(y, x);
}
double Position::dot(const Position p) const {
    return p.x*x + p.y*y;
}

Angle position_to_angle(const Position p) {
    double r = p.norm();
    if (r < DUB_TOL) { r = DUB_TOL;}  // Prevent divide-by-zero
    return Angle(atan2(p.y, p.x), p.y/r, p.x/r);
}

Angle Angle::negate() const {
    return Angle(-a, -sina, cosa);
}  // method (negate)
Angle& Angle::negate_in_place() {
    a    *= -1;
    sina *= -1;
    return *this;
}

Angle Angle::rotate(Angle other) const {
    // Implementation of Angle(c) = Angle(a) + Angle(b):
    // c = a + b
    //
    // [cos(c)] = [cos(b) -sin(b)] * [cos(a)]
    // [sin(c)]   [sin(b)  cos(b)]   [sin(a)]
    return Angle(
        a + other.a,
        sina*other.cosa + cosa*other.sina,
        cosa*other.cosa - sina*other.sina
    );
}  // method (rotate)

Angle& Angle::rotate_in_place(const Angle other) {
    a += other.a;
    double sina_old = sina;
    sina *= other.cosa;
    sina += cosa*other.sina;
    cosa *= other.cosa;
    cosa -= sina_old*other.sina;
    return *this;
}  // method (rotate_in_place)

Turn::Turn() : s(State()), r(0.), nf_id(0) {}
Turn::Turn(State s, double r) : s(s), r(r), nf_id(0) {}
Turn::Turn(State s, double r, unsigned short int nf_id) : s(s), r(r), nf_id(nf_id) {}
std::string Turn::to_string() const {
    return "Turn[" + s.to_string() + " | " + "r=" + std::to_string(r) + " nf_id=" + std::to_string(nf_id) + "]";
}

DubinsPath connect_csc(Turn t0, Turn tf, double tol) {
    Position pc0   = pw2pc(t0.s, t0.r);
    Position pcf   = pw2pc(tf.s, tf.r);
    Position dpc0f = pcf - pc0;
    double dc0f    = dpc0f.norm();

    double sin_tan_csc_num = t0.r - tf.r;

    if (dc0f < tol + fabs(sin_tan_csc_num)) {
        // No path connects these turns
        return DubinsPath(t0, tf, true, std::vector<double>(), State(), State(), DUB_INF);
    }

    double sign0  = copysign(1., t0.r);
    double signf  = copysign(1., tf.r);
    double psi_tan_csc = dpc0f.angle() + asin(sin_tan_csc_num / dc0f);

    // x1/x2 are beginning/end of straight arc
    State x1_csc  = State(pc0, psi_tan_csc);
    x1_csc.p      = pc2pw(x1_csc, t0.r);
    State x2_csc  = State(pcf, x1_csc.psi);
    x2_csc.p      = pc2pw(x2_csc, tf.r);
    Position dp12 = x2_csc.p - x1_csc.p;

    // Construct tha in place
    return DubinsPath(
        t0, tf, true,
        std::vector<double>({
            wrap_dpsi(psi_tan_csc - t0.s.psi.a, sign0, tol),
            dp12.norm(),
            wrap_dpsi(tf.s.psi.a - psi_tan_csc, signf, tol)
        }), x1_csc, x2_csc
    );
}  // method (connect_csc)

DubinsPath connect_ccc(Turn t0, Turn tf, double tol) {
    Position pc0   = pw2pc(t0.s, t0.r);
    Position pcf   = pw2pc(tf.s, tf.r);
    Position dpc0f = pcf - pc0;
    double dc0f    = dpc0f.norm();

    double sin_tan_csc_num = t0.r - tf.r;
    double r0_abs = fabs(t0.r);
    double rf_abs = fabs(tf.r);

    if (dc0f < tol + fabs(sin_tan_csc_num) || dc0f > 2. + r0_abs + rf_abs) {
        // No CCC path connects these turns
        return DubinsPath(t0, tf, false, std::vector<double>(), State(), State(), DUB_INF);
    }

    double sign0    = copysign(1., t0.r);
    double signf    = copysign(1., tf.r);
    double psic0f   = dpc0f.angle();
    double cos_d1   = law_of_cos(1. + rf_abs, dc0f, 1. + r0_abs);
    double cos_d2   = law_of_cos(1. + r0_abs, dc0f, 1. + rf_abs);
    double psi1_ccc = psic0f + sign0 * (acos(cos_d1) + 0.5*pi);
    double psi2_ccc = psic0f - sign0 * (acos(cos_d2) + 0.5*pi);
    State x1_ccc    = State(pc0, psi1_ccc);
    x1_ccc.p        = pc2pw(x1_ccc, t0.r);
    State x2_ccc    = State(pcf, psi2_ccc);
    x2_ccc.p        = pc2pw(x2_ccc, tf.r);

    // Construct tha in place
    return DubinsPath(
        t0, tf, false,
        std::vector<double>({
            wrap_dpsi(psi1_ccc   - t0.s.psi.a, sign0, tol),
            wrap_dpsi(psi2_ccc   - psi1_ccc,  -sign0, tol),
            wrap_dpsi(tf.s.psi.a - psi2_ccc,   signf, tol)
        }), x1_ccc, x2_ccc
    );
}  // function (connect_ccc)

DubinsPath::DubinsPath() : t0(Turn()), tf(Turn()), csc(true), tha(std::vector<double>({0., 0., 0.})), x1(State()), x2(State()), cost(0.) {}
DubinsPath::DubinsPath(Turn t0, Turn tf, bool csc, std::vector<double> tha, State x1, State x2) : t0(t0), tf(tf), csc(csc), tha(tha), x1(x1), x2(x2) {
    double r_middle = csc ? 1. : -copysign(1., t0.r);
    cost            = t0.r*tha[0] + r_middle*tha[1] + tf.r*tha[2];
}
DubinsPath::DubinsPath(Turn t0, Turn tf, bool csc, std::vector<double> tha, State x1, State x2, double cost) : t0(t0), tf(tf), csc(csc), tha(tha), x1(x1), x2(x2), cost(cost) {}

std::string DubinsPath::to_string() const {
    std::string word_str = csc ? "CSC" : "CCC";
    return "DubinsPath[" + word_str + " | " + t0.to_string() + " -> " + tf.to_string() + "]";
}

double DubinsPath::cmp(const DubinsPath& other) const {
    return cost - other.cost;
}
bool DubinsPath::operator>(const DubinsPath& other) const {
    return cmp(other) > 0;
}
bool DubinsPath::operator<(const DubinsPath& other) const {
    return cmp(other) < 0;
}
bool DubinsPath::operator==(const DubinsPath& other) const {
    return cmp(other) == 0;
}

Position pw2pc(State x, double r) {
    return Position(x.p.x - r*x.psi.sina, x.p.y + r*x.psi.cosa);
}

Position pc2pw(State x, double r) {
    return pw2pc(x, -r);
}

// Trigonometry
double law_of_cos(double side_opp, double side_1, double side_2) {
    return (side_1*side_1 + side_2*side_2 - side_opp*side_opp) / (2 * side_1 * side_2);
}

double wrap_dpsi(double dpsi, double direction, double tol) {
    double dpsi_abs = ufmod(direction * dpsi + tol, 2 * pi) - tol;
    return direction * dpsi_abs;
    // if (dpsi_abs < 0.) {
    //     return 0.;
    // } else {
    //     return direction * dpsi_abs;
    // }
}
double wrap_dpsi(double dpsi, double direction) {return wrap_dpsi(dpsi, direction, 0.);}

double wrap_to_pi(double psi) {
    return ufmod(psi + pi, 2*pi) - pi;
}

DubinsPath connect_configuration(State s0, double r0, State sf, double rf) {
    Turn t0 = Turn(s0, r0);  // +/+
    Turn tf = Turn(sf, rf);

    // Calculate 6 candidate paths and pick shortest
    DubinsPath dp_lsl         = connect_csc(t0, tf);
    DubinsPath *shortest_path = &dp_lsl;
    DubinsPath dp_lrl         = connect_ccc(t0, tf);
    if (dp_lrl < *shortest_path) {shortest_path = &dp_lrl;}
    tf.r *= -1.;
    DubinsPath dp_lsr = connect_csc(t0, tf);
    if (dp_lsr < *shortest_path) {shortest_path = &dp_lsr;}
    t0.r *= -1.;
    DubinsPath dp_rsr = connect_csc(t0, tf);
    if (dp_rsr < *shortest_path) {shortest_path = &dp_rsr;}
    DubinsPath dp_rlr = connect_ccc(t0, tf);
    if (dp_rlr < *shortest_path) {shortest_path = &dp_rlr;}
    tf.r = rf;
    DubinsPath dp_rsl = connect_csc(t0, tf);
    if (dp_rsl < *shortest_path) {shortest_path = &dp_rsl;}

    return *shortest_path;
}  // function (connect_configuration)
DubinsPath connect_configuration(State s0, State sf) {return connect_configuration(s0, 1., sf, 1.);}

double plane_cut(const Position p, const Position v) {
    return v.y*p.x - v.x*p.y;
}
double plane_cut(const Position p, const Position p0, const Position v) {
    return v.y*(p.x - p0.x) - v.x*(p.y - p0.y);
}
double plane_cut(const Position p, const LineSegment l) {
    return l.v.y*(p.x - l.p0.x) - l.v.x*(p.y - l.p0.y);
}  // function (plane_cut)

bool trig_in_range_wrapped(const Position cos_sin, const Position lower, const Position upper, const double tol) {
    Position midpoint = 0.5*(lower + upper);
    if (plane_cut(midpoint, lower) > 0) {
        // Midpoint is in the convex out-of-range region (sweep > pi)
        return !(plane_cut(cos_sin, upper) < -tol && plane_cut(cos_sin, lower, -lower) < -tol);
    } else {
        // Midpoint is in the convex in-range region (sweep <= pi)
        return (plane_cut(cos_sin, lower) <= tol && plane_cut(cos_sin, upper, -upper) <= tol);
    }
}  // function (trig_in_range_wrapped)

void Polygon::calculate_neighborhood() {
    // Calculate circle circumscribing the polygon as follows:
    // Cycle through each pair of vertices (Pi, Pj) and choose the pair furthest
    // apart. The center is the mean (Pi+Pj)/2 and the radius is half the distance.
    Position pf = edges.back().p0;
    pf += edges.back().v;

    Position p1 = pf;  // Start/end points of diameter
    Position p2 = pf;
    rc2         = 0.;  // diameter squared
    double d2;         // Candidate diameter squared
    for (std::vector<LineSegment>::iterator it0 = edges.begin(); it0 != edges.end(); it0++) {
        // First, compare with pf b/c this vertex is not in list of edges
        d2 = (pf - it0->p0).norm2();
        if (rc2 < d2) {
            p1  = it0->p0;
            p2  = pf;
            rc2 = d2;
        }

        // Next, cycle though remaining positions
        for (std::vector<LineSegment>::iterator itf = it0 + 1; itf != edges.end(); itf++) {
            d2 = (itf->p0 - it0->p0).norm2();
            if (rc2 < d2) {
                p1  = it0->p0;
                p2  = itf->p0;
                rc2 = d2;
            }
        }  // for (Cycle through pf)
    }  // for (Cycle through p0)

    pc   = (p1 + p2)/2;  // Compute midpoint of the diameter line segment
    rc2 /= 4;            // Convert diameter squared to radius squared
    rc  = sqrt(rc2);
}  // method (calculate_neighborhood)

bool Polygon::position_inside(const Position p) const {
    // Implementation of the improved winding algorithm [2] based on the
    // original winding algorithm in [1]. I clean up the algorithm so that it makes
    // more sense to me. Following Fig. 5 of Ref. [2], the algorithm is essentially:
    // 2w += sign(x_intercept) * (sign(yf) - sign(y0))
    // Where sign(0) = 0, facilitating the "half" additions specified by [1] and [2].
    // By calculating 2w instead of w, integer math may be used.
    //
    // [1] D.G. Alciatore and R. Miranda, "A winding number and point-in-polygon algorithm." Glaxo Virtual Anatomy
    // Project Research Report, Deparment of Mechanical Engineering, Colorado State University (1995).
    // [2] G. N. Kumar and M. Bangi, "An Extension to Winding Number and Point-in-Polygon Algorithm," IFAC 2018,
    // issn. 2405-8963, pp. 548-553, DOI: https://doi.org/10.1016/j.ifacol.2018.05.092.
    if ((p - pc).norm2() > rc2) {
        // Position is outside neighborhood
        return false;
    }  // If (neighborhood check)

    int w = 0;    // Winding number
    Position p0;  // Position of vertex relative to point
    double xint;  // X-intercept location
    double yf;    // y position at end of edge
    for (const LineSegment &edge : edges) {
        p0 = edge.p0 - p;

        // Calculate X-intercept
        yf = p0.y + edge.v.y;
        if (yf*p0.y < 0) {
            // X-intercept is crossed at an intermediate point
            xint = p0.x - p0.y * (edge.v.x / edge.v.y);
        } else if (p0.y == 0) {
            // X-intercept is crossed at the initial point
            xint = p0.x;
        } else if (yf == 0) {
            // X-intercept is crossed at the terminal point
            xint = p0.x + edge.v.x;
        } else {
            // X-intercept is not crossed
            continue;
        }

        // Increment winding number based on direction of intercept
        if (xint < 0) {
            // Negative X-axis is crossed
            // (this is part of [2] but ignored in [1])
            if (yf < 0) {
                w += 1;
            } else if (yf > 0) {
                w -= 1;
            }  // if / else if
            if (p0.y < 0) {
                w -= 1;
            } else if (p0.y > 0) {
                w += 1;
            }  // if / else if
        } else {
            // Positive X-axis is crossed
            if (yf > 0) {
                w += 1;
            } else if (yf < 0) {
                w -= 1;
            }  // if / else if
            if (p0.y > 0) {
                w -= 1;
            } else if (p0.y < 0) {
                w += 1;
            }  // if / else if
        }  // if / else (sign of X intercept)
    }  // for (cycle through edges)

    return w != 0;
}  // method (is_position_inside)

bool Polygon::line_intersect(const LineSegment ls, const double tol) const {
    if (ls.v.norm2() < DUB_TOL) {
        // The line segment is essentially a point
        return position_inside(ls.p0);
    }  // if (ls is a point)

    // Neighborhood check
    double d_proj = (pc - ls.p0).dot(ls.v) / ls.v.norm2();
    if (d_proj > 1) {d_proj = 1;} else if (d_proj < 0) {d_proj = 0;}
    if ((ls.p0-pc + d_proj*ls.v).norm2() > rc2) {
        // Line segment is outside neighborhood
        return false;
    }  // if (neighborhood check)

    // Cycle through edges to check for an intersection
    double tol1 = 1. - tol;
    Position dp;
    double ai_dot_v0;
    double d;   
    for (const LineSegment &edge : edges) {
        dp = edge.p0 - ls.p0;
        ai_dot_v0 = edge.v.x*ls.v.y - edge.v.y*ls.v.x;
        if (fabs(ai_dot_v0) < tol + DUB_TOL) {
            // Lines are parallel -> no intersection exists
            continue;
        }  // if (lines are parellel)

        // Check if line intersection is on the line segment
        d = (edge.v.x*dp.y - edge.v.y*dp.x)/ai_dot_v0;
        if (d < tol || tol1 < d) {
            // Intersection is out-of-bounds for line segment
            continue;
        }

        // Check if line intersection is on the polygon edge
        d = (ls.v.x*dp.y - ls.v.y*dp.x)/ai_dot_v0;
        if (d < tol || tol1 < d) {
            // Intersection is out-of-bounds for polygon edge
            continue;
        }

        // Collision is valid
        return true;
    }  // for (cycle through polygon edges)

    return false;  // No collision with any polygon edge
}  // method (line_intersect)

bool Polygon::turn_intersect(const Position pc0, const double r0, const Position p0, const Position p1, const double tol) const {
    // Find whether a circular turn intersects the polygon. The turn is parameterized by:
    // pc0 -- center of circle
    // r0  -- radius of circle
    // p0  -- initial point of turn (CCW)
    // p1  -- terminal point of turn (CCW)
    double r_sum = fabs(r0) + rc;
    if ((pc0 - pc).norm2() > r_sum*r_sum) {
        // Turn is outside the neighborhood of polygon
        return false;
    }  // if (neighborhood check)

    double r02 = r0*r0;

    Position dp;    // Point on edge relative to turning circle center
    Position pint;  // Intersection position on polygon edge
    double ddp2;    // Magnitude squared of dp
    double v2;      // Magnitude squared of change in position
    double dproj;   // Distance along edge closest to circle center
    double ddproj;  // Distance between dproj and the intersection distance
    for (const LineSegment &edge : edges) {
        // Find the point along the edge closest to the turning circle center
        dp    = edge.p0 - pc0;
        v2    = edge.v.norm2();
        dproj = -dp.dot(edge.v) / v2;
        dp   += dproj*edge.v;
        ddp2  = dp.norm2();
        if (ddp2 > r02) {
            // Edge does not intersect turning circle anywhere
            continue;
        }  // if (check if edge intersects turning circle)

        // Check if the two intersections lie both on the edge and on the turning arc
        ddproj = sqrt((r02 - ddp2)/v2);
        std::vector<double> d_int_vals = {dproj + ddproj, dproj - ddproj};
        for (const double d_int : d_int_vals) {
            if (d_int < 0. || 1. < d_int) {
                // Intersection point is outside polygon edge bounds
                continue;
            }

            pint = edge.p0 + d_int * edge.v;
            bool valid_intersect = r0 < 0 ? trig_in_range_wrapped(pint - pc0, p1 - pc0, p0 - pc0, -tol) : trig_in_range_wrapped(pint - pc0, p0 - pc0, p1 - pc0, -tol);
            if (!valid_intersect) {continue;}

            // Intersection point is valid
            return true;
        }  // for (cycle through two intersections)
    }  // for (cycle through edges)

    return false;  // No collision with any polygon edge
}  // method (turn_intersect)

Polygon Polygon::convex_circumscription() const {
    // Returns a new polygon with a subset of the vertices circumscribing this polygon
    std::vector<LineSegment> edges_convex;
    edges_convex.reserve(edges.size());  // As many or fewer edges

    Position vert_prev = edges.back().p0;
    Position vert_next;
    for (std::vector<LineSegment>::const_iterator it = edges.begin(); it != edges.end(); it++) {
        vert_next = it->p0 + it->v;
        if (plane_cut(it->p0, LineSegment(vert_prev, vert_next - vert_prev)) < 0) {
            // This vertex is inside the convex hull -> skip
            continue;
        }
        edges_convex.push_back(*it);
        vert_prev = it->p0;
    }  // for (cycle through edges)

    return Polygon(edges_convex);
}

bool ConvexHull::point_intersect(const Position p, const double tol) const {
    // If Point is outside the half-plane for any edge, it is not in
    // the convex hull.
    for (const LineSegment &edge : edges) {
        if (plane_cut(p, edge) > 0.) {
            return false;
        }  // if
    }  // for

    return true;
}  // method (point_intersect)

bool ConvexHull::line_intersect(const Position p0, const Position pf, const double tol) const {
    // The line segment is parameterized by its initial and final position
    // p0 and pf. There are no extreme values on the line to minimize the plane cut value,
    // so only the end points are candidates for violating the half-plane of each edge.
    // Therefore, if either point violates the half-plane, the half-plane is violated. If neither
    // point violates one of the half-planes, the line segment falls outside the convex hull.
    for (const LineSegment &edge : edges) {
        if (plane_cut(p0, edge) > 0 && plane_cut(pf, edge) > 0) {
            return false;
        }  // if
    }  // for

    return true;
}  // method (line_intersect)
