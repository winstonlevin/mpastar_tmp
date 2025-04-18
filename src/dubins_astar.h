#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <utility>

#include "dubins.h"
#include "astar.h"
#include "rational_trig.h"

// Configuration classes
struct DangerZone : Circle {
    // Fields
    double cost_gain;

    // Constructors
    DangerZone() : Circle(), cost_gain(0.) {};
    DangerZone(Position pc, double rc, double cost_gain) : Circle(pc, rc), cost_gain(cost_gain) {};

    std::string to_string() const {
        std::stringstream ss;
        ss << "SoftNoFlyZone[(" << pc.x << ", " << pc.y << "), r=" << rc << ", Q=" << cost_gain << "]";
        return ss.str();
    };

    double cost(Position p) const {
        p -= pc;  // Normalized distance offset
        p /= rc;
        double d2 = p.norm2();
        return d2 < 1. ? cost_gain * (1./d2 - 1.) : 0.;
    };
};

struct Connection {
    float cost;  // Total path length of Dubins path (Inf <=> Infeasible)
    bool left0;  // True -> initial left turn, False -> initial right turn
    bool csc;    // True -> CSC path, False -> CCC path
    bool leftf;  // True -> final left turn, False -> final right turn

    Connection() : cost(0.), left0(true), csc(true), leftf(true) {};
    Connection(float cost, bool left0, bool csc, bool leftf) : cost(cost), left0(left0), csc(csc), leftf(leftf) {};

    std::string to_string() const {
        std::stringstream ss;
        ss << "Connection[";
        ss << left0 ? "L" : "R";
        ss << csc ? "S" : "C";
        ss << leftf ? "L" : "R";
        ss << ", cost=" << cost << "]";
        return ss.str();
    };
};

// Graph search classes
struct SampledNFAstar {
    // Fields
    State x0;                           // Initial state constraint
    Position pf;                        // Terminal position constraint
    std::vector<Circle> nf;             // Circular no-fly zones to avoid
    std::vector<Polygon> poly;          // Polygonal obstacles to avoid

    double v, mu_max;                   // Dimensional parameters ([v] = length/time, [mu_max] = 1/time)
    double tol;                         // Tolerance used for generic inequality constraints
    int n_samples;                      // Number of heading samples to use around no-fly zones

    std::vector<AstarNode> nodes;       // Nodes storing turns
    double dpsi;                        // Uniform spacing between heading sampled
    std::vector<Angle> psi_samples;
    std::vector<Turn> turns;

    // Constructors
    SampledNFAstar();
    SampledNFAstar(State x0, Position pf);
    SampledNFAstar(State x0, Position pf, std::vector<Circle> nf);
    SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<Polygon> poly);
    SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max);
    SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<Polygon> poly, double v, double mu_max);
    SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, int n_samples);
    SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<Polygon> poly, int n_samples);
    SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max, int n_samples);
    SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<Polygon> poly, double v, double mu_max, int n_samples);

    // Object methods
    std::string to_string() const;

    bool straight_collision(const State s0, const double d_max, const Circle nf_k) const;
    bool turn_collision(const State s0, const double r0, const double dpsi_max, const Circle nf_k) const;
    DubinsPath connect_edge(const State s0, const double r0, const State sf, const double rf) const;

    void initialize_nodes();
    void sample_boundary(Circle boundary, unsigned short int nf_id, std::vector<Angle> sampled_headings);
    void successor_function(AstarNode *node);
    AstarNode solve_astar();
}; // struct (SampledNFAstar)

unsigned short int project_monotonic_increasing(const double x_val, const std::vector<double> x_grid);
unsigned short int project_uniform_spacing(double x_offset, double x_spacing);

struct DenseTreeNFAstarNodeAux {
    AstarNode n;
    unsigned short int psi_idx;
    State x;

    DenseTreeNFAstarNodeAux();
    DenseTreeNFAstarNodeAux(AstarNode n, unsigned short int psi_idx, State x);

    std::string to_string() const;
};

// struct CollisionCheck {
//     bool collision;
//     bool collision_init;

//     CollisionCheck();
//     CollisionCheck(bool collision);

//     void update_collision(bool collision);

//     std::string to_string() const;
// };

template<typename T>
struct Initializable {
    T value;
    bool init;

    Initializable() : init(false) {};
    Initializable(T value) : value(value), init(true) {};

    void update(T new_value) {
        value = new_value;
        init  = true;
    };

    std::string to_string() const {
        std::stringstream ss;
        ss << "Initializable[";
        ss << (init ? std::to_string(value) : "?");
        ss << "]";
        return ss.str(); 
    };
};

struct DenseHeuristic {
    float heuristic;
    float distance;
    bool terminal;

    DenseHeuristic() : heuristic(0.), distance(0.), terminal(false) {};
    DenseHeuristic(float distance) : heuristic(distance), distance(distance), terminal(false) {};
    DenseHeuristic(float heuristic, float distance) : heuristic(heuristic), distance(distance), terminal(false) {};
    DenseHeuristic(float heuristic, float distance, bool terminal) : heuristic(heuristic), distance(distance), terminal(terminal) {};

    std::string to_string() const {
        return "DenseHeuristic[h=" + std::to_string(heuristic) + ", d=" + std::to_string(distance) + "]";
    };
};

struct DenseTreeNFAstar {
    // Fields
    State x0;                   // Initial state (p nondimensionalized by turn radius)
    Position pf;                // Terminal position (p nondimensionalized by turn radius)
    std::vector<Circle> nf;  // No-fly zones (pc and rc nondimensionalized by turn radius)
    std::vector<DangerZone> dz;  // Danger zones (pc and rc nondimensionalized by turn radius)

    double v, mu_max;                   // Dimensional information ([v] = length/time, [mu_max] = 1/time -> r_min = v/mu_max)
    double tol;                         // Tolerance used for inequality constraints
    double x_min, x_max, y_min, y_max;
    double dpsi, dp;                                             // Spacing of discretized states dpsi, dp=dx=dy
    double dpf_tol;                                              // Tolerance for terminal distance to final position
    unsigned short int n_psi_samples, n_x_samples, n_y_samples;  // Discretization of (x,y,psi)
    unsigned short int x0_idx, y0_idx, psi0_idx, xf_idx, yf_idx; // Discretized idces associated with x0 and pf
    unsigned int next_node_id;                                   // Incrementor for adding unique node IDs

    std::vector<double> x_samples, y_samples, psi_samples, cos_samples, sin_samples;
    std::unordered_map<unsigned int, DenseTreeNFAstarNodeAux> nodes;
    std::vector<std::vector<Initializable<bool>>> collisions;         // 2-D (x, y) array holding whether collision takes place

    // Constructors
    DenseTreeNFAstar();
    DenseTreeNFAstar(State x0, Position pf);
    DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf);
    DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz);
    DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max);
    DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, double v, double mu_max);
    DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, unsigned short int n_psi_samples);
    DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, unsigned short int n_psi_samples);
    DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max, unsigned short int n_psi_samples);
    DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, double v, double mu_max, unsigned short int n_psi_samples);

    // Object methods
    std::string to_string() const;
    // unsigned int idces2id(unsigned short int x_idx, unsigned short int y_idx, unsigned short int psi_idx);
    bool check_collision(unsigned short int x_idx, unsigned short int y_idx);
    void initialize_nodes();
    DenseHeuristic heuristic_fun(const State x) const;
    void add_successor(AstarNode *node, State &x1, unsigned short int psi1_idx);
    void successor_function(AstarNode *node);
    AstarNode solve_astar();
}; // struct (DenseTreeNFAstar)

struct RationalAngle {
    short int dsina, dcosa;
    short int fac;
    double a;
    bool within_tol;

    RationalAngle() : dsina(0), dcosa(1), a(0.), within_tol(true) {};
    RationalAngle(short int dsina, short int dcosa, short int fac) : dsina(dsina), dcosa(dcosa), fac(fac), a(atan2(dsina, dcosa)), within_tol(true) {};
    RationalAngle(short int dsina, short int dcosa, short int fac, double a) : dsina(dsina), dcosa(dcosa), fac(fac), a(a), within_tol(true) {};

    std::string to_string() const {
        std::stringstream ss;
        ss << "RationalAngle";
        ss << (within_tol ? "[" : "(");
        ss << a;
        ss << (within_tol ? "]" : ")");
        return ss.str();
    };  // method (to_string)

    void check_tol(double den, const double tol) {
        den *= den;  // Convert to den^2
        within_tol = abs(dcosa*dcosa + dsina*dsina - den) < tol*den;
    };
};  // struct (RationalAngle)

struct DenseNFAstarNodeAux {
    AstarNode n;
    unsigned short int x_idx, y_idx, psi_idx;

    DenseNFAstarNodeAux() : n(AstarNode()), x_idx(0), y_idx(0), psi_idx(0) {};
    DenseNFAstarNodeAux(AstarNode n, unsigned short int x_idx, unsigned short int y_idx, unsigned short int psi_idx) : n(n), x_idx(x_idx), y_idx(y_idx), psi_idx(psi_idx) {};

    std::string to_string() const {
        std::stringstream ss;
        ss << "DenseNFAstarNodeAux[" << n.to_string() << ", (" << x_idx << ", " << y_idx << ", " << psi_idx << ")]";
        return ss.str();
    };
};

struct Transition {
    std::vector<std::pair<short int, short int>> dxy_idx;  // Series of (dx,dy) values. Intermediate values check for collisions, terminal values are for child
    unsigned short int psi_idx;             // child psi_idx (corresponding to p_idx + dp_idx[0] + ... + dp_idx[N])
    float edge_cost;                        // Cost to transition to this child
    double arc_length;                      // Arc length to transition to this child (equivalent to edge_cost for Dubins)

    Transition() : psi_idx(0), edge_cost(0.f), arc_length(0.) {};
    Transition(unsigned short int psi_idx, double arc_length) : psi_idx(psi_idx), edge_cost(static_cast<float>(arc_length)), arc_length(arc_length) {};

    std::string to_string() const {
        std::stringstream ss;
        ss << "Transition[{(dx,dy)}=";
        for (const std::pair<short int, short int> &xy_pair : dxy_idx) {
            ss << "(" << xy_pair.first << "," << xy_pair.second << "), ";
        }
        ss << "psi=" << psi_idx <<", d=" << edge_cost << "]";
        return ss.str();
    };
};

std::vector<std::vector<Transition>> generate_successor_template(const std::vector<RationalAngle> &psir_samples, const unsigned short int radius, const double k);

struct DenseNFAstar {
    // Fields
    State x0;                    // Initial state (p nondimensionalized by turn radius)
    Position pf;                 // Terminal position (p nondimensionalized by turn radius)
    std::vector<Circle> nf;      // No-fly zones (pc and rc nondimensionalized by turn radius)
    std::vector<DangerZone> dz;  // Danger zones (pc and rc nondimensionalized by turn radius)
    std::vector<Polygon> poly;   // Polygon obstacles (distances nondimensionalized by turn radius)

    double control_cost_gain;           // "k" in the path cost: J = int(1 + k/2 * u^2, 0, tf)
    double v, mu_max;                   // Dimensional information ([v] = length/time, [mu_max] = 1/time -> r_min = v/mu_max)
    double tol;                         // Tolerance used for inequality constraints
    double x_min, x_max, y_min, y_max;  // Bounds for px, py

    // The positions and heading are sampled to result in repeated positions (i.e. position values are rational numbers)
    // The number of possible positions per turn radius (Np/R) is equal to:
    //                Np/R = psi_fidelity * dp_fidelity
    // Where psi_fidelity is the denominatory of the rational trig function
    // (e.g. 5 -> cos(tha), sin(tha) in {-1, -4/5, -3/5, 0, 3/5, 4/5, 1}) [12 samples]
    // And dp_fidelity is the number of turn radii each straight arc traverses
    //
    // The tolerance in the graph search for finding the initial position (since it is a reverse search) is:
    // R * 
    unsigned short int psi_fidelity;    // Denominator value for rational trig (generates heading samples)
    unsigned short int dp_fidelity;     // Multiple of denominator value (1/dp_fidelity is fraction of turn radius applied to straight arcs)
    unsigned short int radius;        // Stores diameter of turning circle
    double dpf_tol;                     // Tolerance for terminal distance to final position
    unsigned short int n_psi_samples, n_x_samples, n_y_samples;  // Discretization of (x,y,psi)
    unsigned short int x0_idx, y0_idx, psi0_idx, xf_idx, yf_idx; // Discretized idces associated with x0 and pf

    std::vector<double> x_samples, y_samples;
    std::vector<RationalAngle> psir_samples;  // Heading samples, but stores r*sin(tha) and r*cos(tha)
    RationalTrig rat_trig;
    std::vector<std::vector<Transition>> successor_template;  // Stores straight/left/right succession as function of heading index
    std::unordered_map<unsigned int, DenseNFAstarNodeAux> nodes;
    std::vector<std::vector<Initializable<std::pair<bool, double>>>> collisions;         // 2-D (x, y) array holding (1) whether collision takes place, (2) danger value

    // Constructors
    DenseNFAstar();
    DenseNFAstar(State x0, Position pf);
    DenseNFAstar(State x0, Position pf, unsigned short int psi_fidelity, unsigned short int dp_fidelity);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, std::vector<Polygon> poly);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, unsigned short int psi_fidelity, unsigned short int dp_fidelity);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, unsigned short int psi_fidelity, unsigned short int dp_fidelity);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, std::vector<Polygon> poly, unsigned short int psi_fidelity, unsigned short int dp_fidelity);
    DenseNFAstar(State x0, Position pf, double v, double mu_max);
    DenseNFAstar(State x0, Position pf, double v, double mu_max, unsigned short int psi_fidelity, unsigned short int dp_fidelity);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, double v, double mu_max);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, std::vector<Polygon> poly, double v, double mu_max);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max, unsigned short int psi_fidelity, unsigned short int dp_fidelity);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, double v, double mu_max, unsigned short int psi_fidelity, unsigned short int dp_fidelity);
    DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, std::vector<Polygon> poly, double v, double mu_max, unsigned short int psi_fidelity, unsigned short int dp_fidelity);

    // Object methods
    std::string to_string() const;
    unsigned int idces2id(const unsigned short int x_idx, const unsigned short int y_idx, const unsigned short int psi_idx) const;
    std::pair<bool, double> check_collision(unsigned short int x_idx, unsigned short int y_idx);
    DenseHeuristic heuristic_fun(const int x1_idx, const int y1_idx, const int psi1_idx) const;
    void initialize_nodes();
    void successor_function(AstarNode *node);
    std::vector<DenseNFAstarNodeAux> construct_path(unsigned int node_id) const;
    AstarNode solve_astar();
}; // struct (DenseNFAstar)