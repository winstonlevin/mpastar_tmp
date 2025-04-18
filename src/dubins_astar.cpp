#include "dubins_astar.h"
#include <iostream>
#include <iomanip>
#include <numbers>

using namespace std::numbers;

const double DUB_INF = std::numeric_limits<double>::infinity();
const float FLT_INF = std::numeric_limits<float>::infinity();
const double DUB_TOL = std::numeric_limits<double>::epsilon();

// Graph search structures
SampledNFAstar::SampledNFAstar() : x0(State()), pf(Position()), v(1.), mu_max(1.), n_samples(12), tol(1E-6) {}
SampledNFAstar::SampledNFAstar(State x0, Position pf) : x0(x0), pf(pf), v(1.), mu_max(1.), n_samples(12), tol(1E-6) {}
SampledNFAstar::SampledNFAstar(State x0, Position pf, std::vector<Circle> nf) : x0(x0), pf(pf), nf(nf), v(1.), mu_max(1.), n_samples(12), tol(1E-6) {}
SampledNFAstar::SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<Polygon> poly) : x0(x0), pf(pf), nf(nf), poly(poly), v(1.), mu_max(1.), n_samples(12), tol(1E-6) {}
SampledNFAstar::SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max) : x0(x0), pf(pf), nf(nf), v(v), mu_max(mu_max), n_samples(12), tol(1E-6) {}
SampledNFAstar::SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<Polygon> poly, double v, double mu_max) : x0(x0), pf(pf), nf(nf), poly(poly), v(v), mu_max(mu_max), n_samples(12), tol(1E-6) {}
SampledNFAstar::SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, int n_samples) : x0(x0), pf(pf), nf(nf), v(1.), mu_max(1.), n_samples(n_samples), tol(1E-6) {}
SampledNFAstar::SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<Polygon> poly, int n_samples) : x0(x0), pf(pf), nf(nf), poly(poly), v(1.), mu_max(1.), n_samples(n_samples), tol(1E-6) {}
SampledNFAstar::SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max, int n_samples) : x0(x0), pf(pf), nf(nf), v(v), mu_max(mu_max), n_samples(n_samples), tol(1E-6) {}
SampledNFAstar::SampledNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<Polygon> poly, double v, double mu_max, int n_samples) : x0(x0), pf(pf), nf(nf), poly(poly), v(v), mu_max(mu_max), n_samples(n_samples), tol(1E-6) {}

std::string SampledNFAstar::to_string() const {
    return "SampledNFAstar[" + x0.to_string() + " -> " + pf.to_string() + " | " + std::to_string(nf.size()) + " NFs ]";
}

bool SampledNFAstar::straight_collision(const State s0, const double d_max, const Circle nf_k) const {
    // Return true if intersects no-fly zone nf_k, false otherwise
    Position dp = nf_k.pc - s0.p;
    double d_proj = dp.x * s0.psi.cosa + dp.y * s0.psi.sina;
    if (d_proj < 0.) {
        d_proj = 0.;
    } else if (d_proj > d_max) {
        d_proj = d_max;
    }  // if/else if (saturate d_proj)

    // Change from pc - p0 to pc - p_proj
    dp.x -= (s0.psi.cosa * d_proj);
    dp.y -= (s0.psi.sina * d_proj);
    return dp.norm2() + tol < nf_k.rc2;
}  // method (straight_collision)

bool SampledNFAstar::turn_collision(const State s0, const double r0, const double dpsi_max, const Circle nf_k) const {
    double sign_turn = copysign(1., r0);

    Position dp      = nf_k.pc - s0.p;  // Distance from center of turn
    dp.x            += (r0*s0.psi.sina);    // to center of no-fly zone
    dp.y            -= (r0*s0.psi.cosa);

    double d2_min = nf_k.rc + fabs(r0);
    if (dp.x*dp.x + dp.y*dp.y + tol > d2_min*d2_min) {
        // Turn does not overlap no-fly zone -> cannot intersect
        return false;
    }

    double dpsi_hat = wrap_dpsi(atan2(sign_turn*dp.x, -sign_turn*dp.y) - s0.psi.a, sign_turn, tol);
    if (sign_turn*dpsi_hat > sign_turn*dpsi_max) {
        dpsi_hat = dpsi_max;
    }

    // Change from pc - p0 to pc - p_proj
    dpsi_hat += s0.psi.a;  // convert dpsi to psi
    dp.x     -= r0 * sin(dpsi_hat);
    dp.y     += r0 * cos(dpsi_hat);

    return dp.x*dp.x + dp.y*dp.y + tol < nf_k.rc2;
}  // method (turn_collision)

DubinsPath SampledNFAstar::connect_edge(const State s0, const double r0, const State sf, const double rf) const {
    // Initialize Priority Queue to avoid unnecessary no-fly zone collisions
    std::priority_queue<DubinsPath, std::vector<DubinsPath>, std::greater<DubinsPath>> dubins_paths;

    // Calculate Dubins paths between configurations
    double r0_rev  = -copysign(1., r0);
    double rf_rev  = -copysign(1., rf);

    dubins_paths.push(connect_csc(Turn(s0, r0),     Turn(sf, rf),     tol));
    dubins_paths.push(connect_csc(Turn(s0, r0),     Turn(sf, rf_rev), tol));
    dubins_paths.push(connect_csc(Turn(s0, r0_rev), Turn(sf, rf),     tol));
    dubins_paths.push(connect_csc(Turn(s0, r0_rev), Turn(sf, rf_rev), tol));

    if (fabs(r0_rev - rf_rev) < 1) {
        // r0/rf and r0_rev/rf_rev share sign
        dubins_paths.push(connect_ccc(Turn(s0, r0),     Turn(sf, rf),     tol));
        dubins_paths.push(connect_ccc(Turn(s0, r0_rev), Turn(sf, rf_rev), tol));
    } else {
        // r0/rf_rev and r0_rev/rf share sign
        dubins_paths.push(connect_ccc(Turn(s0, r0),     Turn(sf, rf_rev),     tol));
        dubins_paths.push(connect_ccc(Turn(s0, r0_rev), Turn(sf, rf),         tol));
    }  // if/else (which CCC arcs are active)

    // Evaluate Dubins in order of optimality, checking no-fly zone violation
    DubinsPath d_path;

    do {
        // Select Dubins Path with lowest cost ------------------------------------- //
        d_path = dubins_paths.top();
        dubins_paths.pop();

        if (std::isinf(d_path.cost)) {
            // No valid path exists, break and return infinite cost path
            break;
        }

        // Ensure no collisions occur ---------------------------------------------- //
        for (const Circle &nf_k : nf) {
            // Check middle turn or straight
            // (Checked first because straight is faster to check)
            if (d_path.csc) {
                if (straight_collision(d_path.x1, d_path.tha[1], nf_k)) {
                    d_path.cost = DUB_INF;
                    continue;
                }
            } else {
                if (turn_collision(d_path.x1, -copysign(1., d_path.t0.r), d_path.tha[1], nf_k)) {
                    d_path.cost = DUB_INF;
                    continue;
                }
            }

            // Check initial turn
            if (turn_collision(d_path.t0.s, d_path.t0.r, d_path.tha[0], nf_k)) {
                d_path.cost = DUB_INF;
                continue;
            }

            // Check terminal turn
            if (turn_collision(d_path.x2, d_path.tf.r, d_path.tha[2], nf_k)) {
                d_path.cost = DUB_INF;
                continue;
            }
        }  // for (iterate through no-fly zones)

        for (const Polygon &polyi : poly) {
            // Check middle turn or straight
            // (Checked first because straight is faster to check)
            if (d_path.csc) {
                if (polyi.line_intersect(LineSegment(d_path.x1.p, d_path.x2.p - d_path.x1.p), tol)) {
                    d_path.cost = DUB_INF;
                    continue;
                }
            } else {
                double r1 = d_path.t0.r < 0 ? -1. : 1.;
                if (polyi.turn_intersect(pw2pc(d_path.x1, r1), r1, d_path.x1.p, d_path.x2.p, tol)) {
                    d_path.cost = DUB_INF;
                    continue;
                }
            }

            // Check initial turn
            if (polyi.turn_intersect(d_path.t0.pc(), d_path.t0.r, d_path.t0.s.p, d_path.x1.p, tol)) {
                d_path.cost = DUB_INF;
                continue;
            }

            // Check terminal turn
            if (polyi.turn_intersect(d_path.tf.pc(), d_path.tf.r, d_path.x2.p, d_path.tf.s.p, tol)) {
                d_path.cost = DUB_INF;
                continue;
            }
        }  // for (iterate through polygons)

        // No collision occured, this path is both fastest AND collision free -> break and accept
        break;
    } while (!dubins_paths.empty());

    return d_path;
}  // method (connect_edge)

void SampledNFAstar::initialize_nodes() {
    // Generate psi, cos(psi), sin(psi) samples for No-Fly zones (evenly spaced, not double counting +-pi)
    psi_samples.reserve(n_samples);

    double psi0  = -pi;
    dpsi = n_samples < 1 ? 2*pi : 2*pi / static_cast<double>(n_samples);

    for (int idx = 0; idx < n_samples; idx++) {
        psi_samples.emplace_back(Angle(psi0));
        psi0 += dpsi;
    }  // for (generate psi samples)

    // CREATE NODES FOR A* GRAPH:
    //
    // Heuristic k  is: cartesian distance from pk to pf
    // Edge cost ij is: Shortest Dubins path from si to sj (with fixed initial/terminal turn sign)
    // Successors to i: Nodes j where Dubins path does not violate no-fly zones
    unsigned int id          = 0;  // counters for unique Node ID and unique circular boundary ID
    unsigned short int nf_id = 0;
    const size_t n_nodes = 2 + 1 + (nf.size() + poly.size()) * n_samples * 2;
    nodes.reserve(n_nodes);
    turns.reserve(n_nodes);

    // Initial node (initial state)
    nodes.emplace_back(AstarNode(id++, static_cast<float>((pf - x0.p).norm()), true, false));
    turns.emplace_back(Turn(x0, 1., nf_id++));

    // Terminal node (point -> rf = 0)
    nodes.emplace_back(AstarNode(id++, std::vector<AstarNode::AstarEdge>(), 0., false, true));
    turns.emplace_back(Turn(State(pf, 0., 0., 1.), 0., nf_id++));

    // ----------------------------------------------- //
    // SAMPLING ALONG CIRCULAR OBSTACLE BOUNDARY:      //
    // sample the states tangent to the circle at the  //
    // sampled heading values                          //
    // ----------------------------------------------- //
    Position p0;
    double rt;
    for (const Circle &nf_k : nf) {
        rt   = nf_k.rc < 1. ? 1. : nf_k.rc;  // Saturate turn radius at 1 (min. turn radius)

        for (const Angle &psi0 : psi_samples) {
            // Sampled position (on circular boundary)
            p0 = Position(nf_k.pc.x + nf_k.rc*psi0.sina, nf_k.pc.y - nf_k.rc*psi0.cosa);
    
            // Include CCW turn
            nodes.emplace_back(AstarNode(id++, static_cast<float>((pf - p0).norm()), false, false));
            turns.emplace_back(Turn(State(p0, psi0), rt, nf_id));
    
            // Include CW turn (same postion -> same heuristic)
            nodes.emplace_back(AstarNode(id++, nodes.back().heuristic, false, false));
            turns.emplace_back(Turn(State(p0, psi0.rotate_180()), -rt, nf_id));
        }  // for (iterate over psi samples)
        nf_id++;
    }  // for (iterate over no-fly zones)

    // ------------------------------------------------- //
    // SAMPLING ALONG POLYGON OBSTACLE BOUNDARY:         //
    // sample the convex vertices at the sampled         //
    // heading values (values may be repeated for non-   //
    // convex polygons)                                  //
    // ------------------------------------------------- //
    Position vert_prev;
    Position vert_next;
    Position v_in;
    Angle psi_min, psi_max;  // Minimum and maximum headings without colliding with polygon edges
    float heur;
    for (const Polygon &polyi : poly) {
        vert_prev = polyi.edges.back().p0;
        v_in      = polyi.edges.back().v;
        for (std::vector<LineSegment>::const_iterator edge_it = polyi.edges.begin(); edge_it != polyi.edges.end(); v_in = edge_it->v, edge_it++) {
            vert_next = edge_it->p0 + edge_it->v;
            if (plane_cut(edge_it->p0, LineSegment(vert_prev, vert_next - vert_prev)) < 0) {
                // This vertex is not on the convex boundary of the polygon
                // [Note: this operation assumes verticies are in CCW order]
                continue;
            }
            vert_prev = edge_it->p0;
            heur      = static_cast<float>((pf - edge_it->p0).norm());

            // Determine valid heading samples ------------------------- //
            Position psi0p;
            for (const Angle &psi0 : psi_samples) {
                psi0p = angle_to_position(psi0);
                if (
                    (plane_cut(psi0p, v_in) <= 0 && 0 <= plane_cut(psi0p, edge_it->v)) // Valid CCW heading
                    || (plane_cut(psi0p, edge_it->v) <= 0 && 0 <= plane_cut(psi0p, v_in))  // Valid CW heading
                ) {
                    nodes.emplace_back(AstarNode(id++, heur, false, false));
                    turns.emplace_back(Turn(State(edge_it->p0, psi0), 1., nf_id));
                }
            }  // for (iterate over possible psi samples)
            nf_id++;
        }  // for (iterate over polygon edges)
    }  // for (iterate over polygon obstacles)
}  // method (initialize_nodes)

void SampledNFAstar::sample_boundary(Circle boundary, unsigned short int nf_id, std::vector<Angle> sampled_headings) {
    unsigned int id = nodes.size();
    double rt   = boundary.rc < 1. ? 1. : boundary.rc;  // Saturate turn radius at 1 (min. turn radius)
    double u_lr[2] = {1., -1.};
    double x0, y0, dx0, dy0;
    float heur;
    for (const Angle &psi0 : sampled_headings) {
        // Calculate heuristic
        x0   = boundary.pc.x + boundary.rc*psi0.sina;
        y0   = boundary.pc.y - boundary.rc*psi0.cosa;
        dx0  = pf.x - x0;
        dy0  = pf.y - y0;
        heur = static_cast<float>(sqrt(dx0*dx0 + dy0*dy0));

        // Include CCW turn
        nodes.emplace_back(AstarNode(id++, heur, false, false));
        turns.emplace_back(Turn(State(x0, y0, psi0), rt, nf_id));

        // Include CW turn
        nodes.emplace_back(AstarNode(id++, heur, false, false));
        turns.emplace_back(Turn(State(x0, y0, psi0.rotate_180()), -rt, nf_id));
    }  // for (iterate over psi samples)
}

void SampledNFAstar::successor_function(AstarNode *node) {
    const Turn *pturn_i = &turns[node->id];

    // ----------------------------------------------------------- //
    //      Connect current STATE to final position                //
    // ----------------------------------------------------------- //
    // Set sign of terminal turn to be same as initial turn
    float edge_cost = static_cast<float>(connect_edge(pturn_i->s, pturn_i->r, turns[1].s, copysign(0., pturn_i->r)).cost);
    if (!std::isinf(edge_cost)) {
        // A valid path to the terminal node exists -> skip intermediate nodes
        node->edges.emplace_back(AstarNode::AstarEdge(&nodes[1], edge_cost));
        return;
    }  // if

    // ----------------------------------------------------------- //
    //      Connect current TURN to subsequent no-fly zones        //
    // ----------------------------------------------------------- //
    std::vector<Turn>::iterator turn_it = turns.begin() + 2;
    for (std::vector<AstarNode>::iterator node_it = nodes.begin() + 2; node_it != nodes.end(); node_it++, turn_it++) {
        if (pturn_i->nf_id == turn_it->nf_id) {
            continue;  // No-fly zones do not self-connect
        }
        edge_cost = static_cast<float>(connect_edge(pturn_i->s, pturn_i->r, turn_it->s, turn_it->r).cost);
        if (!std::isinf(edge_cost)) {
            node->edges.emplace_back(AstarNode::AstarEdge(&*node_it, edge_cost));
        }
    }  // for (inbound node ->j)
}  // method (successor_function)

AstarNode SampledNFAstar::solve_astar() {
    std::priority_queue<AstarNode*, std::vector<AstarNode*>, PointerGEqual> open_list;
    open_list.push(&nodes[0]);  // Initial nodes (1st 2 nodes are open)
    open_list.push(&nodes[1]);
    std::function<void(AstarNode*)> succ_fun = [this](AstarNode *node_i) {return successor_function(node_i);};
    return astar_search(open_list, succ_fun);
}  // method (solve_astar)

// ---------------------------------------------------------------- //
// DENSE (TREE) NO-FLY ZONE                                         //
// ---------------------------------------------------------------- //
unsigned short int project_monotonic_increasing(const double x_val, const std::vector<double> x_grid) {
    unsigned short int idx_min = 0;
    unsigned short int idx_max = static_cast<unsigned short int>(x_grid.size()) - 1;
    unsigned short int idx_new;
    double cost_min            = x_grid[idx_min] - x_val;
    double cost_max            = x_grid[idx_max] - x_val;
    double cost_new;

    // In case x_val is out-of-bounds with min/max value of x_grid
    if (cost_min > 0.) {
        return idx_min;
    } else if (cost_max < 0.) {
        return idx_max;
    }

    // x_val is in-bounds of min/max value of x_grid -> find this value
    // std::cout << "--- idx_min / idx_new / idx_max | cost_min / cost_new / cost_max ---" << std::endl;
    while (idx_max - idx_min > 1) {
        idx_new  = (idx_min + idx_max) / 2;  // Average index
        cost_new = x_grid[idx_new] - x_val;

        // std::cout << "    " << std::setw (7) << idx_min << " / " << std::setw (7) << idx_new << " / " << std::setw (7) << idx_max << " | ";
        // std::cout << std::setw (8) << cost_min << " / " << std::setw (8) << cost_new << " / " << std::setw (8) << cost_max << std::endl;

        if (cost_new > 0.) {
            cost_max = cost_new;
            idx_max = idx_new;
        } else {
            cost_min = cost_new;
            idx_min  = idx_new;
        }
    }

    if (fabs(cost_max) < fabs(cost_min)) {
        return idx_max;
    } else {
        return idx_min;
    }  // if/else (return lower-error index)
}  // function (project_monotonic_increasing)

unsigned short int project_uniform_spacing(double x_offset, double x_spacing) {
    // x_offset:  x_val - x_min
    // x_spacing: uniform change in x across discretization
    //
    // Returns the index i in the discretized array nearest to x_val. The array
    // is assumed to take the form:
    //
    // [x_min + 0*x_spacing, x_min + 1*x_spacing, ..., x_min + (N-1)*x_spacing]
    //
    // Also assumes x_offset < x_max - x_min, or the return index will be
    // out-of-bounds.
    return static_cast<unsigned short int>(round(x_offset/x_spacing));
}

DenseTreeNFAstarNodeAux::DenseTreeNFAstarNodeAux() : psi_idx(0), x(State()) {}
DenseTreeNFAstarNodeAux::DenseTreeNFAstarNodeAux(AstarNode n, unsigned short int psi_idx, State x) : n(n), psi_idx(psi_idx), x(x) {}

std::string DenseTreeNFAstarNodeAux::to_string() const {
    return "DenseTreeNFAstarNodeAux[" + n.to_string() + ", " + std::to_string(psi_idx) + ", " + x.to_string() + "]";
}

// CollisionCheck::CollisionCheck() : collision(false), collision_init(false) {}
// CollisionCheck::CollisionCheck(bool collision) : collision(collision), collision_init(true) {}
// void CollisionCheck::update_collision(bool collision_new) {
//     collision      = collision_new;
//     collision_init = true;
// }  // method (update_collision)
// std::string CollisionCheck::to_string() const {
//     std::string col_str = collision_init ? (collision ? "Yes" : "No") : "?";
//     return "CollisionCheck[" + col_str + "]";
// }  // method (to_string)


DenseTreeNFAstar::DenseTreeNFAstar() : x0(State()), pf(Position()), nf(std::vector<Circle>(0)), v(1.), mu_max(1.), n_psi_samples(12), tol(1E-6), next_node_id(0) {}
DenseTreeNFAstar::DenseTreeNFAstar(State x0, Position pf) : x0(x0), pf(pf), nf(std::vector<Circle>(0)), v(1.), mu_max(1.), n_psi_samples(12), tol(1E-6), next_node_id(0) {}
DenseTreeNFAstar::DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf) : x0(x0), pf(pf), nf(nf), v(1.), mu_max(1.), n_psi_samples(12), tol(1E-6), next_node_id(0) {}
DenseTreeNFAstar::DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max) : x0(x0), pf(pf), nf(nf), v(v), mu_max(mu_max), n_psi_samples(12), tol(1E-6), next_node_id(0) {}
DenseTreeNFAstar::DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, unsigned short int n_psi_samples) : x0(x0), pf(pf), nf(nf), v(1.), mu_max(1.), n_psi_samples(n_psi_samples), tol(1E-6), next_node_id(0) {}
DenseTreeNFAstar::DenseTreeNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max, unsigned short int n_psi_samples) : x0(x0), pf(pf), nf(nf), v(v), mu_max(mu_max), n_psi_samples(n_psi_samples), tol(1E-6), next_node_id(0) {}

std::string DenseTreeNFAstar::to_string() const {
    return "DenseTreeNFAstar[" + x0.to_string() + " -> " + pf.to_string() + " | " + std::to_string(nf.size()) + " NFs ]";
}

bool DenseTreeNFAstar::check_collision(const unsigned short int x_idx, const unsigned short int y_idx) {
    // Cycle through no-fly zones to see if x, y violate a no-fly zone boundary
    Initializable<bool> &collision = collisions[x_idx][y_idx];
    if (collision.init) {
        return collision.value;
    } else {
        double dx, dy;
        const double x = x_samples[x_idx];
        const double y = y_samples[y_idx];

        for (const Circle nf_k : nf) {
            dx = x - nf_k.pc.x;
            dy = y - nf_k.pc.y;
            if (dx*dx + dy*dy < nf_k.rc2) {
                collision.update(true);
                return true;
            }
        }  // for (cycle through no-fly zones)

        collision.update(false);
        return false;
    }  // if (collision_init)
}  // method (check_collision)

DenseHeuristic DenseTreeNFAstar::heuristic_fun(const State s_i) const {
    // The baseline heuristic is the cartesian distance from p -> pf
    Position dp_i = pf - s_i.p;
    DenseHeuristic h = DenseHeuristic(static_cast<float>(dp_i.norm()));

    if (h.distance > 2.) {
        // If dp > 2. -> terminal position is outside of turn radius
        // (no need to calculate turn center)
        return h;
    }

    // Check if terminal state inside turn radius. If so -> replace with PI to account for U-turn.
    // First, convert dp = pf - p0 to dp = pf - pc0 = pf - p0 - dp0c
    if (dp_i.y * s_i.psi.cosa - dp_i.x * s_i.psi.sina < 0.) {
        // Right Turn
        dp_i.x -= s_i.psi.sina;
        dp_i.y += s_i.psi.cosa;
    } else {
        // Left Turn
        dp_i.x += s_i.psi.sina;
        dp_i.y -= s_i.psi.cosa;
    }  // if/else

    // Next, replace with |PI*r|=PI if inside turn radius (dist <= 2, so PI > dist)
    // to account for U-turn required to reach terminal position (admissible)
    if (dp_i.norm2() < 1.) {
        h.heuristic = pi_v<float>;
    } else if (h.distance < dpf_tol) {
        // If within tolerance (and outside turn circle), this state is "terminal"
        h.heuristic = 0.;
        h.terminal  = true;
    }

    return h;
}  // method (heuristic_fun)

void DenseTreeNFAstar::initialize_nodes() {
    // Generate psi, cos(psi), sin(psi) samples (evenly spaced, not double counting +-pi)
    psi_samples.resize(n_psi_samples);
    cos_samples.resize(n_psi_samples);
    sin_samples.resize(n_psi_samples);

    double psi0  = -pi;
    dpsi         = n_psi_samples < 1 ? 0 : 2*pi / static_cast<double>(n_psi_samples);

    for (unsigned short int idx = 0; idx < n_psi_samples; idx++) {
        psi_samples[idx] = psi0;
        cos_samples[idx] = cos(psi0);
        sin_samples[idx] = sin(psi0);
        psi0 += dpsi;
    }  // for (generate psi samples)

    // Generate box bounds for p samples (permit 3 turn radii from outermost points)
    if (pf.x < x0.p.x) {
        x_min = pf.x;
        x_max = x0.p.x;
    } else {
        x_min = x0.p.x;
        x_max = pf.x;
    }
    if (pf.y < x0.p.y) {
        y_min = pf.y;
        y_max = x0.p.y;
    } else {
        y_min = x0.p.y;
        y_max = pf.y;
    }

    double x_min_k, x_max_k, y_min_k, y_max_k;
    for (Circle nf_k : nf) {
        x_min_k = nf_k.pc.x - nf_k.rc;
        x_max_k = nf_k.pc.x + nf_k.rc;
        y_min_k = nf_k.pc.y - nf_k.rc;
        y_max_k = nf_k.pc.y + nf_k.rc;
        if (x_min_k < x_min) {x_min = x_min_k;}
        if (x_max_k > x_max) {x_max = x_max_k;}
        if (y_min_k < y_min) {y_min = y_min_k;}
        if (y_max_k > y_max) {y_max = y_max_k;}
    }  // for (check if no-fly zones change boundary)

    // Generate p samples (spaced dpsi/2 apart within box constraints)
    dp          = dpsi;
    dpf_tol     = dp + tol;

    // add buffer of 3 turn radii (permit U-turns, a.k.a. CCC arcs)
    // (add dpf_tol so that discretization of allowable headings does not create issue)
    double bounds_tol = dpf_tol + 3.;
    x_min -= bounds_tol;
    x_max += bounds_tol;
    y_min -= bounds_tol;
    y_max += bounds_tol;

    n_x_samples = static_cast<unsigned short int>(ceil((x_max - x_min) / dp));
    n_y_samples = static_cast<unsigned short int>(ceil((y_max - y_min) / dp));

    // Fit x_max, y_max to discretization tolerance (for checking bounds)
    x_max = x_min + dp * (n_x_samples - 0.5) - tol;
    y_max = y_min + dp * (n_y_samples - 0.5) - tol;

    x_samples.reserve(n_x_samples);
    double xi = x_min;
    x_samples.emplace_back(xi);
    for (unsigned short int idx = 0; idx < n_x_samples; idx++) {
        xi += dp;
        x_samples.emplace_back(xi);
    }

    y_samples.reserve(n_y_samples);
    xi = y_min;
    for (unsigned short int idx = 0; idx < n_x_samples; idx++) {
        xi += dp;
        y_samples.emplace_back(xi);
    }

    // Project x0, pf onto grid
    x0_idx              = project_uniform_spacing(x0.p.x - x_samples[0], dp);
    y0_idx              = project_uniform_spacing(x0.p.y - y_samples[0], dp);
    double psi0_wrapped = wrap_to_pi(x0.psi.a);
    psi0_idx            = project_uniform_spacing(psi0_wrapped - psi_samples[0], dpsi);
    xf_idx              = project_uniform_spacing(pf.x - x_samples[0], dp);
    yf_idx              = project_uniform_spacing(pf.y - y_samples[0], dp);

    // CREATE NODES FOR A* GRAPH:
    //
    // Heuristic k  is: cartesian distance from pk to pf
    // Edge cost ij is: dpsi
    // Successors to i: terminal states after Left/Straight/Right, projected onto grid
    unsigned int max_nodes = n_x_samples * n_y_samples * n_psi_samples;

    nodes.rehash(max_nodes);      // Ensure Hash map is of sufficient size to store all nodes without re-hashing

    // Reserve 2D collisions array
    collisions.resize(n_x_samples);
    for (std::vector<Initializable<bool>>& collisions_x : collisions) {
        collisions_x.resize(n_y_samples);
    }

    // Insert initial node
    DenseHeuristic h0 = heuristic_fun(x0);
    nodes.emplace(next_node_id, DenseTreeNFAstarNodeAux(AstarNode(next_node_id, h0.heuristic, true, h0.terminal), psi0_idx, x0));                   
    next_node_id++;
}  // method (initialize_nodes)

// unsigned int DenseTreeNFAstar::idces2id(unsigned short int x_idx, unsigned short int y_idx, unsigned short int psi_idx) {
//     return psi_idx + n_psi_samples * y_idx + (n_psi_samples * n_x_samples) * x_idx;
// }

void DenseTreeNFAstar::add_successor(AstarNode *node, State &x1, unsigned short int psi1_idx) {
    if (x1.p.x < x_min || x1.p.x > x_max || x1.p.y < y_min || x1.p.y > y_max) {
        return;  // Position is out of bounds -> skip
    }
    const unsigned short int x1_idx = project_uniform_spacing(x1.p.x - x_samples[0], dp);
    const unsigned short int y1_idx = project_uniform_spacing(x1.p.y - y_samples[0], dp);

    if (!check_collision(x1_idx, y1_idx)) {
        // Valid successor -> add node (if not already added)
        DenseHeuristic h0 = heuristic_fun(x1);
        std::pair<std::unordered_map<unsigned int, DenseTreeNFAstarNodeAux>::iterator, bool> node_return = nodes.emplace(
            next_node_id, DenseTreeNFAstarNodeAux(AstarNode(next_node_id, h0.heuristic, false, h0.terminal), psi1_idx, x1)
        );
        next_node_id++;
        node->edges.emplace_back(AstarNode::AstarEdge(&node_return.first->second.n, static_cast<float>(dpsi)));
    }
}

void DenseTreeNFAstar::successor_function(AstarNode *node) {
    const DenseTreeNFAstarNodeAux &aux_info = nodes[node->id];
    State x1 = aux_info.x;
    unsigned short int psi1_idx;

    node->edges.reserve(3);

    // Straight path
    x1.p.x += dpsi * x1.psi.cosa;
    x1.p.y += dpsi * x1.psi.sina;
    psi1_idx = aux_info.psi_idx;
    add_successor(node, x1, psi1_idx);

    // Left turn path
    x1 = aux_info.x;
    psi1_idx = aux_info.psi_idx == n_psi_samples - 1 ? 0 : aux_info.psi_idx + 1;  // Wrap heading
    x1.psi.a = psi_samples[psi1_idx];      // modify heading
    x1.psi.cosa = cos_samples[psi1_idx];
    x1.psi.sina = sin_samples[psi1_idx];
    x1.p.x += x1.psi.sina - aux_info.x.psi.sina;  //  modify position with turn
    x1.p.y += aux_info.x.psi.cosa - x1.psi.cosa;
    add_successor(node, x1, psi1_idx);

    // Right turn path
    x1 = aux_info.x;
    psi1_idx = aux_info.psi_idx == 0 ? n_psi_samples - 1 : aux_info.psi_idx - 1;  // Wrap heading
    x1.psi.a = psi_samples[psi1_idx];      // modify heading
    x1.psi.cosa = cos_samples[psi1_idx];
    x1.psi.sina = sin_samples[psi1_idx];
    x1.p.x -= x1.psi.sina - aux_info.x.psi.sina;  //  modify position with turn
    x1.p.y -= aux_info.x.psi.cosa - x1.psi.cosa;
    add_successor(node, x1, psi1_idx);
}

AstarNode DenseTreeNFAstar::solve_astar() {
    std::priority_queue<AstarNode*, std::vector<AstarNode*>, PointerGEqual> open_list;
    open_list.push(&nodes[0].n);  // Initial node (corresponding to initial state)
    std::function<void(AstarNode*)> succ_fun = [this](AstarNode *node_i) {return successor_function(node_i);};
    return astar_search(open_list, succ_fun);
}  // method (solve_astar)

// ---------------------------------------------------------------- //
// DENSE (GRID) NO-FLY ZONE                                         //
// ---------------------------------------------------------------- //
std::vector<std::vector<Transition>> generate_successor_template(const std::vector<RationalAngle> &psir_samples, const unsigned short int radius, const double k) {
    std::vector<std::vector<Transition>> transitions = std::vector<std::vector<Transition>>();

    // Calculate transitions based on Rational psir sampling
    const unsigned short int psi_idx_max = static_cast<unsigned short int>(psir_samples.size() - 1);
    const unsigned short int n_transitions = static_cast<unsigned short int>(psi_idx_max + 1);
    unsigned short int psi0_idx, psi1_idx;

    // double r_dubins       = static_cast<unsigned short int>(round(sqrt(k/2)));
    double turn_gain_full = 1. + k/2;

    transitions.resize(n_transitions);
    short int dx_prev, dx, ddx, dx_total;
    short int dy_prev, dy, ddy, dy_total;
    unsigned short int arc_length_total;
    Transition *tr_ptr;

    for (unsigned short int psi_idx = 0; psi_idx < n_transitions; psi_idx++) {
        transitions.reserve(3);

        // Straight (go forward radius/fac integer amount, checking intermediate
        // values to no-fly zone violations)
        dx_total         = psir_samples[psi_idx].dcosa / psir_samples[psi_idx].fac;
        dy_total         = psir_samples[psi_idx].dsina / psir_samples[psi_idx].fac;
        arc_length_total = radius / psir_samples[psi_idx].fac;

        transitions[psi_idx].emplace_back(Transition(psi_idx, arc_length_total));

        dx_prev = 0;
        dy_prev = 0;
        transitions[psi_idx].back().dxy_idx.reserve(arc_length_total);
        for (unsigned short int arc_length = 1; arc_length <= arc_length_total; arc_length++) {
            dx = static_cast<short int>(round(static_cast<double>(arc_length)/static_cast<double>(arc_length_total) * static_cast<double>(dx_total)));
            dy = static_cast<short int>(round(static_cast<double>(arc_length)/static_cast<double>(arc_length_total) * static_cast<double>(dy_total)));

            ddx = dx - dx_prev;
            ddy = dy - dy_prev;

            if (ddx == 0 && ddy == 0) {continue;}

            transitions[psi_idx].back().dxy_idx.emplace_back(std::pair<short int, short int>(ddx, ddy));

            // Advance initial values
            dx_prev = dx;
            dy_prev = dy;
        }

        // Left turn (psi_idx INCREASING until psi1_idx is reached where heading is within tolerance)
        transitions[psi_idx].emplace_back(Transition(
            psi_idx,
            0.
        ));

        psi0_idx          = psi_idx;
        psi1_idx          = psi0_idx == psi_idx_max ? 0 : psi0_idx + 1;
        do {
            ddx                = psir_samples[psi1_idx].dsina - psir_samples[psi0_idx].dsina;
            ddy                = psir_samples[psi0_idx].dcosa - psir_samples[psi1_idx].dcosa;

            transitions[psi_idx].back().dxy_idx.emplace_back(std::pair<short int, short int>(ddx, ddy));

            // Advance initial values
            psi0_idx = psi1_idx;
            psi1_idx = psi0_idx == psi_idx_max ? 0 : psi0_idx + 1;
        } while(!psir_samples[psi0_idx].within_tol);
        tr_ptr = &transitions[psi_idx].back();
        tr_ptr->psi_idx = psi0_idx;
        tr_ptr->arc_length = radius*fabs(
            wrap_to_pi(psir_samples[psi0_idx].a - psir_samples[psi_idx].a)
        );
        tr_ptr->edge_cost = static_cast<float>(tr_ptr->arc_length*turn_gain_full);

        // Right turn (psi_idx DECREASING until psi1_idx is reached where heading is within tolerance)
        transitions[psi_idx].emplace_back(Transition(
            psi_idx,
            0.
        ));

        psi0_idx = psi_idx;
        psi1_idx = psi0_idx == 0 ? psi_idx_max : psi0_idx - 1;
        do {
            ddx                = psir_samples[psi0_idx].dsina - psir_samples[psi1_idx].dsina;
            ddy                = psir_samples[psi1_idx].dcosa - psir_samples[psi0_idx].dcosa;

            transitions[psi_idx].back().dxy_idx.emplace_back(std::pair<short int, short int>(ddx, ddy));

            // Advance initial values
            psi0_idx = psi1_idx;
            psi1_idx = psi0_idx == 0 ? psi_idx_max : psi0_idx - 1;
        } while(!psir_samples[psi0_idx].within_tol);
        tr_ptr = &transitions[psi_idx].back();
        tr_ptr->psi_idx = psi0_idx;
        tr_ptr->arc_length = radius*fabs(
            wrap_to_pi(psir_samples[psi0_idx].a - psir_samples[psi_idx].a)
        );
        tr_ptr->edge_cost = static_cast<float>(tr_ptr->arc_length*turn_gain_full);
    }  // for (each psi transition)

    return transitions;
}  // function (generate_successor_template)

const double TOL_DEFAULT                      = 1E-6;
const unsigned short int PSI_FIDELITY_DEFAULT = 5;
const unsigned short int DP_FIDELITY_DEFAULT  = 2;

DenseNFAstar::DenseNFAstar() : x0(State()), pf(Position()), control_cost_gain(0.), v(1.), mu_max(1.), psi_fidelity(PSI_FIDELITY_DEFAULT), dp_fidelity(DP_FIDELITY_DEFAULT), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf) : x0(x0), pf(pf), control_cost_gain(0.), v(1.), mu_max(1.), psi_fidelity(PSI_FIDELITY_DEFAULT), dp_fidelity(DP_FIDELITY_DEFAULT), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, unsigned short int psi_fidelity, unsigned short int dp_fidelity) : x0(x0), pf(pf), control_cost_gain(0.), v(1.), mu_max(1.), psi_fidelity(psi_fidelity), dp_fidelity(dp_fidelity), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf) : x0(x0), pf(pf), nf(nf), control_cost_gain(0.), v(1.), mu_max(1.), psi_fidelity(PSI_FIDELITY_DEFAULT), dp_fidelity(DP_FIDELITY_DEFAULT), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz) : x0(x0), pf(pf), nf(nf), dz(dz), control_cost_gain(0.), v(1.), mu_max(1.), psi_fidelity(PSI_FIDELITY_DEFAULT), dp_fidelity(DP_FIDELITY_DEFAULT), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, std::vector<Polygon> poly) : x0(x0), pf(pf), nf(nf), dz(dz), poly(poly), control_cost_gain(0.), v(1.), mu_max(1.), psi_fidelity(PSI_FIDELITY_DEFAULT), dp_fidelity(DP_FIDELITY_DEFAULT), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, double v, double mu_max) : x0(x0), pf(pf), control_cost_gain(0.), v(v), mu_max(mu_max), psi_fidelity(psi_fidelity), dp_fidelity(dp_fidelity), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, double v, double mu_max, unsigned short int psi_fidelity, unsigned short int dp_fidelity) : x0(x0), pf(pf), control_cost_gain(0.), v(v), mu_max(mu_max), psi_fidelity(psi_fidelity), dp_fidelity(dp_fidelity), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max) : x0(x0), pf(pf), nf(nf), control_cost_gain(0.), v(v), mu_max(mu_max), psi_fidelity(PSI_FIDELITY_DEFAULT), dp_fidelity(DP_FIDELITY_DEFAULT), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, double v, double mu_max) : x0(x0), pf(pf), nf(nf), dz(dz), control_cost_gain(0.), v(v), mu_max(mu_max), psi_fidelity(PSI_FIDELITY_DEFAULT), dp_fidelity(DP_FIDELITY_DEFAULT), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, std::vector<Polygon> poly, double v, double mu_max) : x0(x0), pf(pf), nf(nf), dz(dz), poly(poly), control_cost_gain(0.), v(v), mu_max(mu_max), psi_fidelity(PSI_FIDELITY_DEFAULT), dp_fidelity(DP_FIDELITY_DEFAULT), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, unsigned short int psi_fidelity, unsigned short int dp_fidelity) : x0(x0), pf(pf), nf(nf), control_cost_gain(0.), v(1.), mu_max(1.), psi_fidelity(psi_fidelity), dp_fidelity(dp_fidelity), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, unsigned short int psi_fidelity, unsigned short int dp_fidelity) : x0(x0), pf(pf), nf(nf), dz(dz), control_cost_gain(0.), v(1.), mu_max(1.), psi_fidelity(psi_fidelity), dp_fidelity(dp_fidelity), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, std::vector<Polygon> poly, unsigned short int psi_fidelity, unsigned short int dp_fidelity) : x0(x0), pf(pf), nf(nf), dz(dz), poly(poly), control_cost_gain(0.), v(1.), mu_max(1.), psi_fidelity(psi_fidelity), dp_fidelity(dp_fidelity), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, double v, double mu_max, unsigned short int psi_fidelity, unsigned short int dp_fidelity) : x0(x0), pf(pf), nf(nf), control_cost_gain(0.), v(v), mu_max(mu_max), psi_fidelity(psi_fidelity), dp_fidelity(dp_fidelity), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, double v, double mu_max, unsigned short int psi_fidelity, unsigned short int dp_fidelity) : x0(x0), pf(pf), nf(nf), dz(dz), control_cost_gain(0.), v(v), mu_max(mu_max), psi_fidelity(psi_fidelity), dp_fidelity(dp_fidelity), tol(TOL_DEFAULT) {}

DenseNFAstar::DenseNFAstar(State x0, Position pf, std::vector<Circle> nf, std::vector<DangerZone> dz, std::vector<Polygon> poly, double v, double mu_max, unsigned short int psi_fidelity, unsigned short int dp_fidelity) : x0(x0), pf(pf), nf(nf), dz(dz), poly(poly), control_cost_gain(0.), v(v), mu_max(mu_max), psi_fidelity(psi_fidelity), dp_fidelity(dp_fidelity), tol(TOL_DEFAULT) {}

std::string DenseNFAstar::to_string() const {
    std::stringstream ss;
    ss << "DenseNFAstar[(" << x0.p.x << ", " << x0.p.y << ", " << x0.psi.a << ") -> ";
    ss << "(" << pf.x << ", " << pf.y << ") | ";
    ss << nf.size() << " NFs, " << dz.size() << " DZs" << "]";
    return ss.str();
}  // method (DenseNFAstar::to_string)

std::pair<bool, double> DenseNFAstar::check_collision(const unsigned short int x_idx, const unsigned short int y_idx) {
    // Cycle through no-fly zones to see if x, y violate a no-fly zone boundary
    Initializable<std::pair<bool, double>> &collision = collisions[x_idx][y_idx];

    if (!collision.init) {
        // Check no-fly zone collisions
        bool collided = false;
        const Position p = Position(x_samples[x_idx], y_samples[y_idx]);
        for (const Circle &nf_k : nf) {
            if ((p - nf_k.pc).norm2() < nf_k.rc2) {
                collided = true;
                break;
            }
        }  // for (cycle through circular obstacles)

        if (!collided) {
            for (const Polygon &poly_i : poly) {
                if (poly_i.position_inside(p)) {
                    collided = true;
                    break;
                }
            }
        }  // for (cycle through polygon obstacles)

        // Check danger zone values
        double danger = 0.;
        for (const DangerZone dz_k : dz) {
            danger += dz_k.cost(p);
        }

        collision.update(std::pair<bool, double>({collided, danger}));
    }  // if (collision_init)

    return collision.value;
}  // method (check_collision)

unsigned int DenseNFAstar::idces2id(const unsigned short int x_idx, const unsigned short int y_idx, const unsigned short int psi_idx) const {
    return psi_idx + n_psi_samples * y_idx + n_psi_samples * n_y_samples * x_idx;
}

void DenseNFAstar::initialize_nodes() {
    rat_trig = RationalTrig(psi_fidelity*dp_fidelity);
    unsigned short int n_rat_trig = static_cast<unsigned short int>(rat_trig.ns.size());

    // Set tolerances / spacing based on trig approximation
    radius      = dp_fidelity*psi_fidelity;
    double dp   = 1./radius;
    dpf_tol     = psi_fidelity + tol;

    // Number of psi samples given by rational trig (but rat trig only includes Quadrant I, i.e. [0, 90] degrees -> 4 angles per ratio (except 0/90 deg)):
    // NOTE: RationalTrig includes (ns, nc) pairs in order of ascending heading from 0 to 90, and includes 0 as its first value and 90 as its last
    n_psi_samples = static_cast<unsigned short int>(4*rat_trig.ns.size()) - 2;
    psir_samples.reserve(n_psi_samples);
    RationalAngle *psir_back_it;
    std::vector<double> psir_proj_double = std::vector<double>();  // For projecting psi0 onto attainable heading
    std::vector<unsigned short int> psir_proj_idces = std::vector<unsigned short int>();
    psir_proj_double.reserve(n_psi_samples);
    psir_proj_idces.reserve(n_psi_samples);

    // Quadrant III angles [include -pi, -pi/2]
    for (unsigned short int idx = 0; idx < n_rat_trig; idx++) {
        psir_samples.emplace_back(RationalAngle(-rat_trig.ns[idx], -rat_trig.nc[idx], rat_trig.fac[idx], rat_trig.tha[idx] - pi));
        psir_back_it = &psir_samples.back();
        psir_back_it->check_tol(radius, tol);
        if (psir_back_it->within_tol) {
            psir_proj_double.emplace_back(psir_back_it->a);
            psir_proj_idces.emplace_back(psir_samples.size() - 1);
        }
    }

    // Quadrant IV angles [exclude -pi/2, -0]
    for (unsigned short int idx = n_rat_trig - 2; idx > 0; idx--) {
    // for (unsigned short int idx = 1; idx < n_rat_trig - 1; idx++) { // <- reverse
        psir_samples.emplace_back(RationalAngle(-rat_trig.ns[idx], rat_trig.nc[idx], rat_trig.fac[idx], -rat_trig.tha[idx]));
        psir_back_it = &psir_samples.back();
        psir_back_it->check_tol(radius, tol);
        if (psir_back_it->within_tol) {
            psir_proj_double.emplace_back(psir_back_it->a);
            psir_proj_idces.emplace_back(psir_samples.size() - 1);
        }
    }

    // Quadrant I angles [include 0, pi/2]
    for (unsigned short int idx = 0; idx < n_rat_trig; idx++) {
        psir_samples.emplace_back(RationalAngle(rat_trig.ns[idx], rat_trig.nc[idx], rat_trig.fac[idx], rat_trig.tha[idx]));
        psir_back_it = &psir_samples.back();
        psir_back_it->check_tol(radius, tol);
        if (psir_back_it->within_tol) {
            psir_proj_double.emplace_back(psir_back_it->a);
            psir_proj_idces.emplace_back(psir_samples.size() - 1);
        }
    }

    // Quadrant II angles [exclude pi/2, pi]
    for (unsigned short int idx = n_rat_trig - 2; idx > 0; idx--) {
    // for (unsigned short int idx = 1; idx < n_rat_trig - 1; idx++) { // <- reverse
        psir_samples.emplace_back(RationalAngle(rat_trig.ns[idx], -rat_trig.nc[idx], rat_trig.fac[idx], pi - rat_trig.tha[idx]));
        psir_back_it = &psir_samples.back();
        psir_back_it->check_tol(radius, tol);
        if (psir_back_it->within_tol) {
            psir_proj_double.emplace_back(psir_back_it->a);
            psir_proj_idces.emplace_back(psir_samples.size() - 1);
        }
    }

    // Generate box bounds for p samples (permit 3 turn radii from outermost points)
    if (pf.x < x0.p.x) {
        x_min = pf.x;
        x_max = x0.p.x;
    } else {
        x_min = x0.p.x;
        x_max = pf.x;
    }
    if (pf.y < x0.p.y) {
        y_min = pf.y;
        y_max = x0.p.y;
    } else {
        y_min = x0.p.y;
        y_max = pf.y;
    }

    double x_min_k, x_max_k, y_min_k, y_max_k;
    for (Circle nf_k : nf) {
        x_min_k = nf_k.pc.x - nf_k.rc;
        x_max_k = nf_k.pc.x + nf_k.rc;
        y_min_k = nf_k.pc.y - nf_k.rc;
        y_max_k = nf_k.pc.y + nf_k.rc;
        if (x_min_k < x_min) {x_min = x_min_k;}
        if (x_max_k > x_max) {x_max = x_max_k;}
        if (y_min_k < y_min) {y_min = y_min_k;}
        if (y_max_k > y_max) {y_max = y_max_k;}
    }  // for (check if no-fly zones change boundary)

    // add buffer of 3 turn radii (permit U-turns, a.k.a. CCC arcs)
    double bounds_tol = 3.;
    x_min -= bounds_tol;
    x_max += bounds_tol;
    y_min -= bounds_tol;
    y_max += bounds_tol;

    // Generate p samples (spaced 1 integer unit apart within box constraints)
    n_x_samples = static_cast<unsigned short int>(ceil((x_max - x_min)*radius));
    n_y_samples = static_cast<unsigned short int>(ceil((y_max - y_min)*radius));

    // Fit x_max, y_max to discretization tolerance (for checking bounds)
    x_max = x_min + (n_x_samples - 0.5)/radius - tol;
    y_max = y_min + (n_y_samples - 0.5)/radius - tol;

    x_samples.reserve(n_x_samples);
    double xi = x_min;
    x_samples.emplace_back(xi);
    for (unsigned short int idx = 0; idx < n_x_samples; idx++) {
        xi += dp;
        x_samples.emplace_back(xi);
    }

    y_samples.reserve(n_y_samples);
    xi = y_min;
    for (unsigned short int idx = 0; idx < n_y_samples; idx++) {
        xi += dp;
        y_samples.emplace_back(xi);
    }

    // Snap initial state / terminal position to grid
    x0_idx              = project_uniform_spacing(x0.p.x - x_samples[0], dp);
    y0_idx              = project_uniform_spacing(x0.p.y - y_samples[0], dp);
    double psi0_wrapped = wrap_to_pi(x0.psi.a);
    psi0_idx            = psir_proj_idces[project_monotonic_increasing(psi0_wrapped, psir_proj_double)];
    xf_idx              = project_uniform_spacing(pf.x - x_samples[0], dp);
    yf_idx              = project_uniform_spacing(pf.y - y_samples[0], dp);

    // Derive template for succession
    successor_template = generate_successor_template(psir_samples, radius, control_cost_gain);

    // CREATE NODES FOR A* GRAPH:
    //
    // Heuristic k  is: cartesian distance from pk to pf
    // Edge cost ij is: dpsi
    // Successors to i: terminal states after Left/Straight/Right, projected onto grid
    unsigned int max_nodes = n_x_samples * n_y_samples * static_cast<unsigned int>(psir_proj_double.size());

    nodes.rehash(max_nodes);      // Ensure Hash map is of sufficient size to store all nodes without re-hashing

    // Reserve 2D collisions array
    collisions.resize(n_x_samples);
    for (std::vector<Initializable<std::pair<bool, double>>>& collisions_x : collisions) {
        collisions_x.resize(n_y_samples);
    }

    // Insert initial node
    DenseHeuristic h0 = heuristic_fun(x0_idx, y0_idx, psi0_idx);
    unsigned int id0  = idces2id(x0_idx, y0_idx, psi0_idx);
    nodes.emplace(id0, DenseNFAstarNodeAux(AstarNode(id0, h0.heuristic, true, h0.terminal), x0_idx, y0_idx, psi0_idx));
}  // method (initialize_nodes)

DenseHeuristic DenseNFAstar::heuristic_fun(const int x1_idx, const int y1_idx, const int psi1_idx) const {
    // The baseline heuristic is the cartesian distance from p -> pf
    double dx = radius*(pf.x - x_samples[x1_idx]);
    double dy = radius*(pf.y - y_samples[y1_idx]);
    DenseHeuristic h = DenseHeuristic(static_cast<float>(sqrt(dx*dx + dy*dy)));

    if (h.distance > 2*radius) {
        // If dp > diameter -> terminal position is outside of turn radius
        // (no need to calculate turn center)
        return h;
    }

    // Check if terminal state inside turn radius. If so -> replace with PI to account for U-turn.
    // First, convert dp = pf - p0 to dp = pf - pc0 = pf - p0 - dp0c
    RationalAngle psir = psir_samples[psi1_idx];

    if (dx * psir.dcosa + dy * psir.dsina < 0.) {
        // True -> target is in BEHIND vehicle
        // If target is BEHIND, U-turn is required so heuristic must be at least PI*R.
        // Additionally, this state is NOT terminal, even if target is close.
        float cost_min = static_cast<float>(pi*radius);
        if (h.heuristic < cost_min) {
            h.heuristic = cost_min;
        }
        return h;
    }

    if (dy * psir.dcosa - dx * psir.dsina < 0.) {
        // Right Turn to target (for CSC path)
        dx -= psir.dsina;
        dy += psir.dcosa;
    } else {
        // Left Turn to target (for CSC path)
        dx += psir.dsina;
        dy -= psir.dcosa;
    }  // if/else

    if (dx*dx + dy*dy < radius*radius) {
        // Target is inside our turn radius. Increase heauristic from dist to |PI*r|.
        // (NOTE: dist <= 2*R, so PI*R > dist)
        // This accounts for U-turn required to reach terminal position (admissible)
        h.heuristic = static_cast<float>(pi*radius);
    } else if (h.distance < dpf_tol) {
        // Within tolerance and target is NOT in our turn radius NOR behind us,
        // so this state is "terminal"
        h.heuristic = 0.;
        h.terminal  = true;
    }

    return h;
}  // method (heuristic_fun)

void DenseNFAstar::successor_function(AstarNode *node) {
    const DenseNFAstarNodeAux &aux_info = nodes[node->id];
    std::pair<std::unordered_map<unsigned int, DenseNFAstarNodeAux>::iterator, bool> node_return;
    unsigned int next_node_id;

    int x1_idx;  // signed to permit subtraction from dx_idx
    int y1_idx;

    node->edges.reserve(3);
    std::vector<unsigned short int>::iterator x1_it, y1_it;
    for (Transition & trans : successor_template[aux_info.psi_idx]) {
        // Follow indices until the child is reached
        x1_idx     = aux_info.x_idx;
        y1_idx     = aux_info.y_idx;
        
        double danger = 0.;
        std::pair<bool, double> collision;
        for (const std::pair<short int, short int> &dxy : trans.dxy_idx) {
            x1_idx += dxy.first;
            y1_idx += dxy.second;

            // Box bounds check
            if (x1_idx < 0 || y1_idx < 0 || x1_idx >= n_x_samples || y1_idx >= n_y_samples) {
                goto continue_successor_loop;  // Out of bounds, skip this node
            }

            collision = check_collision(x1_idx, y1_idx);
            if (collision.first) {
                goto continue_successor_loop;  // Collision occurs, skip this node
            } else {
                danger += collision.second;
            }
        }  // for (series of shifts)

        // Try adding successor (see if it is already added)
        next_node_id = idces2id(x1_idx, y1_idx, trans.psi_idx);
        node_return = nodes.insert({
            next_node_id, DenseNFAstarNodeAux(AstarNode(next_node_id, 0., false, false), x1_idx, y1_idx, trans.psi_idx)
        });
        if (node_return.second) {
            // This node has not been visited before, calculate heuristic
            DenseHeuristic h1 = heuristic_fun(x1_idx, y1_idx, trans.psi_idx);
            node_return.first->second.n.heuristic = h1.heuristic;
            node_return.first->second.n.terminal  = h1.terminal;
        }

        // Add edge
        node->edges.emplace_back(AstarNode::AstarEdge(&node_return.first->second.n, trans.edge_cost + (danger / static_cast<double>(trans.dxy_idx.size()))*trans.arc_length));

        continue_successor_loop:;
    }  // for (successor loop)
}  // method (successor_function)

AstarNode DenseNFAstar::solve_astar() {
    std::priority_queue<AstarNode*, std::vector<AstarNode*>, PointerGEqual> open_list;
    open_list.push(&nodes.begin()->second.n);  // Initial node (corresponding to initial state)
    std::function<void(AstarNode*)> succ_fun = [this](AstarNode *node_i) {return successor_function(node_i);};
    return astar_search(open_list, succ_fun);
}  // method (solve_astar)

std::vector<DenseNFAstarNodeAux> DenseNFAstar::construct_path(unsigned int node_id) const {
    std::vector<DenseNFAstarNodeAux> path = std::vector<DenseNFAstarNodeAux>();
    path.push_back(nodes.at(node_id));
    while (path.back().n.parent) {
        node_id = path.back().n.parent->id;
        path.push_back(nodes.at(node_id));
    }
    return path;
}  // method (construct_path)
