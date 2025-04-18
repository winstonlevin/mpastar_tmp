#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <numbers>

using namespace std::numbers;

bool dub_eq(double a, double b);

double law_of_cos(double side_opp, double side_1, double side_2);
double wrap_dpsi(double dpsi, double direction, double tol);
double wrap_dpsi(double dpsi, double direction);
double wrap_to_pi(double psi);

// Configuration structures
struct Position{
    double x, y;

    Position() : x(0.), y(0.) {};
    // Position(const Position &p) : x(p.x), y(p.y) {};
    Position(double x, double y) : x(x), y(y) {};

    // Operators
    Position operator+(const Position& other) const;
    Position operator+(const double& other) const;
    friend Position operator+(double other, const Position &p);
    Position& operator+=(const Position& other);
    Position& operator+=(const double& other);
    Position operator-(const Position& other) const;
    Position operator-(const double& other) const;
    Position operator-() const;
    friend Position operator-(double other, const Position &p);
    Position& operator-=(const Position& other);
    Position& operator-=(const double& other);
    Position operator*(const double& other) const;
    friend Position operator*(double other, const Position &p);
    Position& operator*=(const double& other);
    Position operator/(const double& other) const;
    friend Position operator/(double other, const Position &p);
    Position& operator/=(const double& other);

    // Methods
    std::string to_string() const;
    double norm2() const;
    double norm() const;
    double angle() const;
    double dot(const Position p) const;
};

struct Angle{
    double a;  // Angle value [rad]
    double sina;  // Sine value (sin(a))
    double cosa;  // Cosine value (cos(a))

    Angle() : a(0.), sina(0.), cosa(1.) {};
    Angle(double a) : a(a), sina(sin(a)), cosa(cos(a)) {};
    Angle(double a, double sina, double cosa) : a(a), sina(sina), cosa(cosa) {};

    std::string to_string() {return "Angle[" + std::to_string(a) + "]";};

    // Trig ratios
    double tana() const {return sina/cosa;};  // Tangent
    double cota() const {return cosa/sina;};  // Cotangent
    double seca() const {return 1/cosa;};     // Secant
    double csca() const {return 1/sina;};     // Cosecant

    // Simple rotations
    Angle negate() const;
    Angle& negate_in_place();
    Angle rotate(Angle other) const;
    Angle& rotate_in_place(const Angle other);
    inline Angle rotate_90() const {return rotate(Angle(0.5*pi, 1., 0.));};
    inline Angle rotate_180() const {return rotate(Angle(pi, 0., -1.));};
    inline Angle rotate_270() const {return rotate(Angle(-0.5*pi, -1., 0.));};
    inline Angle rotate_90_in_place() {return rotate_in_place(Angle(0.5*pi, 1., 0.));};
    inline Angle rotate_180_in_place() {return rotate_in_place(Angle(pi, 0., -1.));};
    inline Angle rotate_270_in_place() {return rotate_in_place(Angle(-0.5*pi, -1., 0.));};

    // Operators
    Angle operator+(const Angle& other) const {return rotate(other);};
    friend Angle operator+(double other, const Angle &psi) {return psi.rotate(other);};
    Angle& operator+=(const Angle& other) {return rotate_in_place(other);};
    Angle operator-() const {return negate();};
    Angle operator-(const Angle& other) const {return rotate(other.negate());};
    friend Angle operator-(double other, const Angle &psi) {return psi.negate().rotate(other);};
    Angle& operator-=(const Angle& other) {return rotate_in_place(other.negate());};

    inline bool in_range_wrapped(const Angle lower, const Angle upper) const {return wrap_to_pi(a - lower.a) <= wrap_to_pi(upper.a - lower.a);};
    inline bool in_range_wrapped(const double lower, const double upper) const {return wrap_to_pi(a - lower) <= wrap_to_pi(upper - lower);};
};

// Chosen equivalence with position:
// Associate Angle[angle, sin(angle), cos(angle)] with Position[cos(angle), sin(angle)]
inline Position angle_to_position(const Angle psi) {return Position(psi.cosa, psi.sina);};
Angle position_to_angle(const Position p);

struct State{
    Position p;
    Angle psi;

    // Constructors
    State() : p(Position()), psi(Angle()) {};

    State(Position p, Angle psi) : p(p), psi(psi) {};

    State(Position p, double psi) : p(p), psi(Angle(psi)) {};
    State(Position p, double psi, double sin_psi, double cos_psi) : p(p), psi(Angle(psi, sin_psi, cos_psi)) {};

    State(double x, double y, Angle psi) : p(Position(x, y)), psi(psi) {};
    State(double x, double y, double psi) : p(Position(x, y)), psi(Angle(psi)) {};
    State(double x, double y, double psi, double sin_psi, double cos_psi) : p(Position(x, y)), psi(Angle(psi, sin_psi, cos_psi)) {};

    // Methods
    std::string to_string() const {
        std::stringstream ss;
        ss << "State[" << p.x << ", " << p.y << ", " << psi.a << "]";
        return ss.str();
    };
};

// Operations on configurations
Position pw2pc(State x, double r);
Position pc2pw(State x, double r);

struct Turn {
    State s;
    double r;
    unsigned short int nf_id;

    Turn();
    Turn(State s, double r);
    Turn(State s, double r, unsigned short int nf_id);

    std::string to_string() const;
    Position pc() {return pw2pc(s, r);};
};

struct DubinsPath {
    Turn t0, tf;             // Save initial turn state and terminal turn state
    State x1, x2;            // State after initial turn/before final turn
    bool csc;                // True -> CSC, False -> CCC
    std::vector<double> tha; // change in heading for each turn, arclength for straight
    double cost;             // Stores total arclength

    DubinsPath();
    DubinsPath(Turn t0, Turn tf, bool csc, std::vector<double> tha, State x1, State x2);
    DubinsPath(Turn t0, Turn tf, bool csc, std::vector<double> tha, State x1, State x2, double cost);

    double cmp(const DubinsPath& other) const;       // Comparison operators compare cost
    bool operator>(const DubinsPath& other) const;   // (for managing priority queue)
    bool operator<(const DubinsPath& other) const;
    bool operator==(const DubinsPath& other) const;

    std::string to_string() const;
};

DubinsPath connect_configuration(State s0, double r0, State sf, double rf);
DubinsPath connect_configuration(State s0, State sf);

DubinsPath connect_csc(Turn t0, Turn tf, double tol);
inline DubinsPath connect_csc(Turn t0, Turn tf) {return connect_csc(t0, tf, 0.);};
DubinsPath connect_ccc(Turn t0, Turn tf, double tol);
inline DubinsPath connect_ccc(Turn t0, Turn tf) {return connect_ccc(t0, tf, 0.);};


struct LineSegment{
    // -------------------------------------- //
    // The Line segment is parameterized as:  //
    // p(d) = p0 + d * v                      //
    // -------------------------------------- //
    Position p0;  // Initial position
    Position v;   // Change in position (pf - p0) of line segment

    LineSegment() : p0(Position()), v(Position()) {};
    LineSegment(Position p0, Position v) : p0(p0), v(v) {};

    std::string to_string() {
        std::stringstream ss;
        ss << "LineSegment[";
        ss << "(" << p0.x << ", " << p0.y << ")";
        ss << " + d(" << v.x << ", " << v.y << ")";
        ss << " | 0 <= d <= 1]";
        return ss.str();
    }
};

struct Circle {
    // Fields
    Position pc;
    double rc, rc2;

    // Constructors
    Circle() : pc(Position()), rc(1.), rc2(1.) {};
    Circle(Position pc) : pc(pc), rc(1.), rc2(1.) {};
    Circle(Position pc, double rc) : pc(pc), rc(rc), rc2(rc*rc) {};

    // Object methods
    std::string to_string() const {
        std::stringstream ss;
        ss << "Circle(" << pc.x << ", " << pc.y << " | r=" << rc << ")";
        return ss.str();
    };
};  // struct (Circle)

struct Polygon{
    // Edges defining the polygon
    std::vector<LineSegment> edges;

    // Center, radius, and radius squared of circle circumscribing the polygon
    Position pc;
    double rc, rc2;

    Polygon() : pc(Position()), rc(0.), rc2(0.) {};
    Polygon(std::vector<LineSegment> edges) : edges(edges) {calculate_neighborhood();};

    std::string to_string() const {
        std::stringstream ss;
        ss << "Polygon[";
        ss << "(" << pc.x << ", " << pc.y << "), r=" << rc;
        ss << " | " << edges.size() << " Edges]";
        return ss.str();
    };

    void calculate_neighborhood();

    bool position_inside(const Position p) const;

    bool line_intersect(const LineSegment ls, const double tol) const;
    inline bool line_intersect(const LineSegment ls) const {return line_intersect(ls, 0.);};

    bool turn_intersect(const Position pc0, const double r0, const Position p0, const Position p1, const double tol) const;
    inline bool turn_intersect(const Position pc0, const double r0, const Position p0, const Position p1) const {return turn_intersect(pc0, r0, p0, p1, 0.);};
    inline bool turn_intersect(const Position pc0, const Position p0, const Position p1, const double tol) const {return turn_intersect(pc0, 0., p0, p1, tol);};
    inline bool turn_intersect(const Position pc0, const Position p0, const Position p1) const {return turn_intersect(pc0, 1., p0, p1, 0.);};

    Polygon convex_circumscription() const;
};

struct ConvexHull{
    std::vector<LineSegment> edges;

    ConvexHull() {};
    ConvexHull(std::vector<LineSegment> edges) : edges(edges) {};

    std::string to_string() const {
        std::stringstream ss;
        ss << "ConvexHull[";
        ss << edges.size();
        ss << " Edges]";
        return ss.str();
    };

    // Collision methods
    bool point_intersect(const Position p, const double tol) const;
    inline bool point_intersect(const Position p) const {return point_intersect(p, 0.);};

    bool line_intersect(const Position p0, const Position pf, const double tol) const;
    inline bool line_intersect(const Position p0, const Position pf) const {return line_intersect(p0, pf, 0.);};
};  // struct (ConvexHull)

// Collision algorithms
double plane_cut(const Position p, const Position v);
double plane_cut(const Position p, const Position p0, const Position v);
double plane_cut(const Position p, const LineSegment l);
bool trig_in_range_wrapped(const Position cos_sin, const Position lower, const Position upper, const double tol);
inline bool trig_in_range_wrapped(const Position cos_sin, const Position lower, const Position upper) {return trig_in_range_wrapped(cos_sin, lower, upper, 0.);};
inline bool trig_in_range_wrapped(const Angle cos_sin, const Position lower, const Position upper, const double tol) {return trig_in_range_wrapped(angle_to_position(cos_sin), lower, upper, tol);};
inline bool trig_in_range_wrapped(const Angle cos_sin, const Position lower, const Position upper) {return trig_in_range_wrapped(cos_sin, lower, upper, 0.);};
