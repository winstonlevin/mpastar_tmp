#include <cmath>
#include <vector>
#include <string>
#include <queue>
#include <functional>

struct PointerLEqual {
    template<typename T>
    bool operator()(T *a, T *b) {
        return (*a) < (*b);
    }
};

struct PointerGEqual {
    template<typename T>
    bool operator()(T *a, T *b) {
        return (*a) > (*b);
    }
};

struct PointerEqual {
    template<typename T>
    bool operator()(T *a, T *b) {
        return (*a) == (*b);
    }
};

// Configuration structures
struct IntegerPosition{
    int x, y;

    // Constructors
    IntegerPosition();
    IntegerPosition(int x, int y);

    // Methods
    std::string to_string();
};

struct Node{
    int i, j, g, h, f, ip, jp;
    bool open, terminal, obstructed;

    // Constructors
    Node();
    Node(int i, int j);
    Node(int i, int j, int g, int h);
    Node(int i, int j, int g, int h, int f);
    Node(int i, int j, int g, int h, int f, bool open);
    Node(int i, int j, int g, int h, int f, bool open, bool terminal);

    // Operators
    int cmp(const Node& other) const;
    bool operator>(const Node& other) const;
    bool operator<(const Node& other) const;
    bool operator==(const Node& other) const;
    
    // Methods
    void update_cost_to_come(int g_new);
    std::string to_string();
};

// Graph searches
struct ManhattanAstar {
    // Fields
    IntegerPosition p0, pf;
    std::vector<IntegerPosition> obstacles;
    int height, width, max_iter, status;
    std::vector<std::vector<Node>> nodes;

    // Constructors
    ManhattanAstar();
    ManhattanAstar(IntegerPosition p0, IntegerPosition pf);
    ManhattanAstar(IntegerPosition p0, IntegerPosition pf, int height, int width);
    ManhattanAstar(IntegerPosition p0, IntegerPosition pf, int height, int width, int max_iter);

    ManhattanAstar(IntegerPosition p0, IntegerPosition pf, std::vector<IntegerPosition> obstacles);
    ManhattanAstar(IntegerPosition p0, IntegerPosition pf, std::vector<IntegerPosition> obstacles, int height, int width);
    ManhattanAstar(IntegerPosition p0, IntegerPosition pf, std::vector<IntegerPosition> obstacles, int height, int width, int max_iter);

    // Helper methods
    void autoset_bounds();

    // Object Methods
    void solve_graph();
    std::vector<Node> construct_path(int i, int j);
    std::string to_string();
};
// std::vector<IntegerPosition> manhattan_astar(IntegerPosition p0, IntegerPosition pf, int width, int height);

// Generic A* Algorithm
struct AstarNode {

    struct AstarEdge {
        AstarNode *child; // pointer to child node
        float edge_cost;  // edge cost to child node
    
        AstarEdge();
        AstarEdge(AstarNode *child, float edge_cost);
    
        std::string to_string();
    };

    unsigned int id;               // Unique ID associated with Node
    AstarNode *parent;             // Current parent's id [No parent -> same as own id]
    std::vector<AstarEdge> edges;  // Edges from this node to children nodes
    float cost_to_come;            // Current cost to come from initial node through parent node
    float heuristic;               // Admissible cost estimate to goal node
    float cost;                    // Stored sum cost_to_come + heuristic
    bool open;                     // Whether node is ``open'' (i.e. needs to be examined by A*)
    bool terminal;                 // Whether node is in goal set (i.e. terminate A* if this node is reached)
    bool req_init;                 // Where node requires initializing edges

    AstarNode();
    AstarNode(unsigned int id, float heuristic, bool open, bool terminal);
    AstarNode(unsigned int id, std::vector<AstarEdge> edges, float heuristic, bool open, bool terminal);

    void init_cost_to_come();

    // Operators
    float cmp(const AstarNode& other) const;        // Comparison operators compare cost
    bool operator>(const AstarNode& other) const;   // (for managing priority queue)
    bool operator<(const AstarNode& other) const;
    bool operator==(const AstarNode& other) const;
    
    // Methods
    bool check_parent(AstarNode *possible_parent, float edge_cost);
    std::string to_string() const;
};

std::vector<AstarNode> construct_path(AstarNode node);
AstarNode astar_search(std::priority_queue<AstarNode*, std::vector<AstarNode*>, PointerGEqual> &open_list, std::function<void(AstarNode*)> &successor_func);
