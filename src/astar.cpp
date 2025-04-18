#include "astar.h"
#include <iostream>

const float FLT_INF = std::numeric_limits<float>::infinity();

// Configuration structures
IntegerPosition::IntegerPosition() {x=0, y=0;}
IntegerPosition::IntegerPosition(int x, int y) : x(x), y(y) {}

std::string IntegerPosition::to_string() {
    return "(" + std::to_string(x) + ", " + std::to_string(y) + ")";
}

Node::Node() : i(0), j(0), g(0), h(0), f(0), ip(-1), jp(-1), open(false), terminal(false), obstructed(false) {}
Node::Node(int i, int j) : i(i), j(j), g(0), h(0), f(0), ip(-1), jp(-1), open(false), terminal(false), obstructed(false) {}
Node::Node(int i, int j, int g, int h) : i(i), j(j), g(g), h(h), f(g + h), ip(-1), jp(-1), open(false), terminal(false), obstructed(false) {}
Node::Node(int i, int j, int g, int h, int f) : i(i), j(j), g(g), h(h), f(f), ip(-1), jp(-1), open(false), terminal(false), obstructed(false) {}
Node::Node(int i, int j, int g, int h, int f, bool open) : i(i), j(j), g(g), h(h), f(f), ip(-1), jp(-1), open(open), terminal(false), obstructed(false) {}
Node::Node(int i, int j, int g, int h, int f, bool open, bool terminal) : i(i), j(j), g(g), h(h), f(f), ip(-1), jp(-1), open(open), terminal(terminal), obstructed(false) {}

// Boolean operators compare estimated cost (g >= g*, h <= h*, f = g + h)
int Node::cmp(const Node& other) const {
    return f - other.f;
}
bool Node::operator>(const Node& other) const {
    return cmp(other) > 0;
}
bool Node::operator<(const Node& other) const {
    return cmp(other) < 0;
}
bool Node::operator==(const Node& other) const {
    return cmp(other) == 0;
}

void Node::update_cost_to_come(int g_new) {
    g = g_new;
    f = g + h;
}

std::string Node::to_string() {
    if (open) {
        return "Node((" + std::to_string(i) + ", " + std::to_string(j) + ") | f=" + std::to_string(f) + ")";
    } else {
        return "Node[(" + std::to_string(i) + ", " + std::to_string(j) + ") | f=" + std::to_string(f) + "]";
    }
}

// Manhattan distance A* implementation, assuming distances are normalized to 1
// and the bottom-left-most position is (0, 0)
void ManhattanAstar::autoset_bounds() {
    // Set bounds to maximum input value with a buffer of 1 (index 0 -> + 2 instead of + 1)
    std::vector<IntegerPosition> points = obstacles;
    points.push_back(p0);
    points.push_back(pf);

    width = 0;
    height = 0;

    int ii, ji;
    for (int idx_point = 0; idx_point < points.size(); idx_point++) {
        ii = points[idx_point].x + 2;
        ji = points[idx_point].y + 2;
        if (height < ii) {
            height = ii;
        }
        if (width < ji) {
            width = ji;
        }
    }
}  // method (autoset_bounds)
ManhattanAstar::ManhattanAstar() : p0(IntegerPosition()), pf(IntegerPosition()), obstacles(std::vector<IntegerPosition>(0)), width(0), height(0), max_iter(0) {}

ManhattanAstar::ManhattanAstar(IntegerPosition p0, IntegerPosition pf) : p0(p0), pf(pf), obstacles(std::vector<IntegerPosition>(0)) {
    autoset_bounds();  // initializes width, height
    max_iter = height*width;
} // constructor

ManhattanAstar::ManhattanAstar(IntegerPosition p0, IntegerPosition pf, int height, int width) : p0(p0), pf(pf), obstacles(std::vector<IntegerPosition>(0)), height(height), width(width), max_iter(height*width) {}

ManhattanAstar::ManhattanAstar(IntegerPosition p0, IntegerPosition pf, int height, int width, int max_iter) : p0(p0), pf(pf), obstacles(std::vector<IntegerPosition>(0)), height(height), width(width), max_iter(max_iter) {}

ManhattanAstar::ManhattanAstar(IntegerPosition p0, IntegerPosition pf, std::vector<IntegerPosition> obstacles) : p0(p0), pf(pf), obstacles(obstacles) {
    autoset_bounds();  // Initializes width, height
    max_iter = width*height;
}

ManhattanAstar::ManhattanAstar(IntegerPosition p0, IntegerPosition pf, std::vector<IntegerPosition> obstacles, int height, int width) : p0(p0), pf(pf), obstacles(obstacles), height(height), width(width), max_iter(width*height) {}

ManhattanAstar::ManhattanAstar(IntegerPosition p0, IntegerPosition pf, std::vector<IntegerPosition> obstacles, int height, int width, int max_iter) : p0(p0), pf(pf), obstacles(obstacles), height(height), width(width), max_iter(max_iter) {}

void ManhattanAstar::solve_graph() {
    // Initialize height x width grid of nodes
    nodes.resize(height);

    for (int idx = 0; idx < height; idx++) {
        nodes[idx].resize(width);
        for (int jdx = 0; jdx < width; jdx++) {
            // Heuristic is distance to terminal idx, rounded up
            int dx = pf.x - idx;
            int dy = pf.y - jdx;

            // Nodes are equivalent to their (idx, jdx) position, all start closed.
            nodes[idx][jdx] = Node(idx, jdx, -1, static_cast<int>(ceil(sqrt(static_cast<float>(dx*dx + dy*dy)))), -1);
        } // for (fill in columns of Nodes)
    } // for (fill in rows of Nodes)

    nodes[p0.x][p0.y].update_cost_to_come(0); // No cost to arrive at initial node
    nodes[pf.x][pf.y].terminal = true;        // Put terminal node in goal set

    // Loop through obstacles to set obstructed==true in grid
    IntegerPosition ob;
    while(!obstacles.empty()) {
        ob = obstacles.back();
        if (ob.x < 0 || ob.x >= height || ob.y < 0 || ob.y >= width) {
            // Obstacle is out of bounds -> does not affect nodes
            continue;
        }
        nodes[ob.x][ob.y].obstructed = true;
        obstacles.pop_back();
    }

    // Implementation of A* algorithm to solve for cost-to-come until a parent to the terminal node is found
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> open_list;

    nodes[p0.x][p0.y].open = true;  // Seed search with initial node
    open_list.push(nodes[p0.x][p0.y]);
    Node node;
    
    // Successors for Manhattan search: +1i, -1i, +1j, -1j
    std::vector<IntegerPosition> successors;
    int x_new, y_new, g_new;
    const int n_successors = 4;
    successors.resize(n_successors);
    int dx[n_successors] = {1, -1, 0,  0};
    int dy[n_successors] = {0,  0, 1, -1};
    
    status = -3;
    for (int iteration = 0; iteration < max_iter; iteration++) {
        // Body of A* loop
        if (open_list.empty()) {
            // No more nodes exist to search
            if (nodes[pf.x][pf.y].g == -1) {
                // No feasible path was found -> graph does not contain any solutions
                status = -1;
            } else {
                // The feasible solution is the optimal solution
                status = 1;
            }
            break;
        } // if (check if open list is empty)

        // Check next node
        node = open_list.top();
        open_list.pop();
        if (node.open) {
            node.open = false;
        } else {
            // This node is already checked but still in queue (b/c cost-to-come was updated)
            continue;
        }  // if/else (handle whether or not node is open)

        if (node.terminal) {
            // Node is in the goal set, A* is complete
            status = 1;
            break;
        }  // if

        // Node is not in goal set, apply successor operator
        std::cout << "i=" << std::to_string(iteration) << ": " << node.to_string() << " | " << open_list.size() << " open nodes" << std::endl;

        std::cout << "Successors:" << std::endl;
        for (int successor_idx = 0; successor_idx < n_successors; successor_idx++) {
            x_new = node.i + dx[successor_idx];
            if (x_new < 0 || x_new >= height) {
                continue;  // x out of bounds
            }  // if
            y_new = node.j + dy[successor_idx];
            if (y_new < 0 || y_new >= width) {
                continue;  // y out of bounds
            }  // if
            if (nodes[x_new][y_new].obstructed) {
                continue;  // (x,y) inside obstacle
            }

            // New position is valid -> check if cost-to-come is reduced (edge cost == 1)
            g_new = node.g + 1;
            std::cout << nodes[x_new][y_new].to_string() << " | g (old) = " << std::to_string(nodes[x_new][y_new].g) << ", g_new = " << std::to_string(g_new) << std::endl;
            if (nodes[x_new][y_new].g == -1 || nodes[x_new][y_new].g > g_new) {
                // Cost-to-come is reduced -> make node parent of successor and add success to open list
                std::cout << "ACCEPT" << std::endl;
                nodes[x_new][y_new].ip = node.i;
                nodes[x_new][y_new].jp = node.j;
                nodes[x_new][y_new].update_cost_to_come(g_new);
                nodes[x_new][y_new].open = true;
                open_list.push(nodes[x_new][y_new]);
            } else{
                std::cout << "REJECT" << std::endl;
            }  // if/else (whether to set node as parent of successor)
        } // for (successors of current node)
    } // for (A* search loop)

    // Print remaining items in open list
    while (!open_list.empty()) {
        node = open_list.top();
        std::cout << node.to_string() << std::endl;
        open_list.pop();
    }

    if (status == -3 && nodes[pf.x][pf.y].g != -1) {
        // Max iterations exceeded, but feasible path is found
        status = 3;
    }
} // method (solve_graph)

std::vector<Node> ManhattanAstar::construct_path(int i, int j) {
    // Constructs path backward from nodes[i][j] through parents and puts them into a vector
    std::vector<Node> path;

    // Build path from terminal node back to beginning
    int ip;
    while (i != -1 && j != -1) {
        std::cout << nodes[i][j].to_string() << std::endl;
        path.push_back(nodes[i][j]);
        ip = nodes[i][j].ip;
        j  = nodes[i][j].jp;
        i  = ip;
    }

    return path;
} // method (construct_path)

std::string ManhattanAstar::to_string() {
    return "ManhattanAstar[" + std::to_string(height) + "x" + std::to_string(width) + " | " + p0.to_string() + "->" + pf.to_string() + "]";
} // method (to_string)


// Generic A* Algorithm
AstarNode::AstarEdge::AstarEdge() : child(nullptr), edge_cost(0.) {}
AstarNode::AstarEdge::AstarEdge(AstarNode *child, float edge_cost) : child(child), edge_cost(edge_cost) {}
std::string AstarNode::AstarEdge::to_string() {
    return "AstarEdge(->" + child->to_string() + " | e=" + std::to_string(edge_cost) + ")";
}

AstarNode::AstarNode() : id(0), parent(nullptr), cost_to_come(FLT_INF), heuristic(FLT_INF), cost(FLT_INF), open(false), terminal(false), req_init(true) {}
AstarNode::AstarNode(unsigned int id, float heuristic, bool open, bool terminal) : id(id), parent(nullptr), heuristic(heuristic), open(open), terminal(terminal), req_init(true) {
    init_cost_to_come();
}
AstarNode::AstarNode(unsigned int id, std::vector<AstarEdge> edges, float heuristic, bool open, bool terminal) : id(id), parent(nullptr), edges(edges), heuristic(heuristic), open(open), terminal(terminal), req_init(false) {
    init_cost_to_come();
}
void AstarNode::init_cost_to_come() {
    if (open) {
        cost_to_come = 0;
        cost         = heuristic;
    } else {
        cost_to_come = FLT_INF;
        cost         = FLT_INF;
    }
}

float AstarNode::cmp(const AstarNode& other) const {
    return cost - other.cost;
}
bool AstarNode::operator>(const AstarNode& other) const {
    return cmp(other) > 0;
}
bool AstarNode::operator<(const AstarNode& other) const {
    return cmp(other) < 0;
}
bool AstarNode::operator==(const AstarNode& other) const {
    return cmp(other) == 0;
}

bool AstarNode::check_parent(AstarNode *possible_parent, float edge_cost) {
    // Resassign parent if it results in lower cost-to-come,
    // return whether parent is reassigned
    float cost_to_come_new = possible_parent->cost_to_come + edge_cost;
    if (cost_to_come_new < cost_to_come) {
        // Replace parent with new parent
        parent       = possible_parent;
        cost_to_come = cost_to_come_new;
        cost         = cost_to_come_new + heuristic;
        open         = true;
        return true;
    } else {
        return false;
    }
}

std::string AstarNode::to_string() const {
    std::string parent_string = parent ? std::to_string(parent->id) : "";
    std::string children_string = req_init ? "?" : "{" + std::to_string(edges.size()) + " children}";
    std::string beg_str;
    std::string end_str;
    if (open) {
        beg_str = "A*Node(";
        end_str = ")";
    } else {
        beg_str = "A*Node[";
        end_str = "]";
    }
    return beg_str + parent_string + "->" + std::to_string(id) + "->" + children_string + " | g=" + std::to_string(cost_to_come) + ", h=" + std::to_string(heuristic) + end_str;
}

std::vector<AstarNode> construct_path(AstarNode node) {
    // Constructs path backward from leaf node back to original parent (the Adam node).
    // The Adam node does not have a parent
    std::vector<AstarNode> path = {node};
    while (node.parent) {
        node = *node.parent;
        path.emplace_back(node);
    }  // while (loop until Adam node is reached)
   return path; 
}  // function (construct_path)

AstarNode astar_search(std::priority_queue<AstarNode*, std::vector<AstarNode*>, PointerGEqual> &open_list, std::function<void(AstarNode*)> &successor_func) {
    // The original A* algorithm, sourced from:
    //
    // Hart, P. E., Nilsson, N. J., and Raphael, B., "A Formal Basis for the Heuristic
    // Determination of Minimum Cost Paths", IEEE Transactions of Systems Science and
    // Cybernetics, Vol. SSC-4, No. 2, 1968.
    // URL: https://doi.org/10.1109/TSSC.1968.300136.

    // Assuming open list is already initialized
    // (A* Step 1)
    std::cout << "Running A*..." << std::endl;

    AstarNode* current_node;

    while (!open_list.empty()) {
        // Check open node with lowest cost
        // (A* Step 2)
        current_node = open_list.top();
        open_list.pop();
        std::cout << "Current Node: " << current_node->to_string() << "(" << std::to_string(open_list.size()) << " more open nodes)" << std::endl;

        if (current_node->open) {
            // Close current node (A* algorithm step 3a/4a)
            current_node->open = false;
        } else {
            // This node is already checked but still in queue (b/c cost-to-come was updated)
            std::cout << "(Already checked node, skipping)" << std::endl;
            continue;
        }  // if/else (conditions for open nodes)

        if (current_node->terminal) {
            // Current node is in goal set, terminal algorithm (A* step 3b)
            std::cout << "(Terminal node found, terminating)" << std::endl;
            break;
        }  // if (terminal node)

        if (current_node->req_init) {
            // Edges outbound from current node are not yet calculated, apply successor function
            successor_func(current_node);
            current_node->req_init = false;
        }

        for (const AstarNode::AstarEdge &edge : current_node->edges) {
            // Apply successor function (A* algorithm step 4b)
            if (edge.child->check_parent(current_node, edge.edge_cost)) {
                open_list.push(edge.child);
            }
        }  // for (iterate through edges)
    }  // while (A* search loop)

    // Return goal node
    return *current_node;
}  // function (astar_search)
