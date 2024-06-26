/*
 * file: graph.h
 * Implementing a graph using adjacency list representation
 * System: Mac using CLion
 * Author: Himanshu Dongre
 */

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <stack>
#include <queue>
#include <deque>
#include <set>
#include <map>
#include <unordered_map>
#include <functional>

using namespace std;

//
// class graph
//
// onne edge between two vertices
// an edge is uniquely defined by its starting
// and ending vertices only
//
template<typename vT, typename wT>
class graph {
private:
    //
    // Adjacency list implementation
    //
    unordered_map<vT, map<vT, wT>> Adj_lst;

    //
    // forms_cycle
    //
    // Private member helper function
    // Using in add_edge to check whether an added edge
    // forms a cycle
    //
    bool forms_cycle(const vT &from, const vT &to, deque<vT> cycle, set<vT> discovered);

public:
    //
    // default constructor
    //
    graph() = default;

    //
    // adding vertices for a graph from a set
    //
    graph(const set<vT> &vertices);

    //
    // operator==
    //
    friend bool operator==(const graph<vT, wT> &G1, const graph<vT, wT> &G2)
    {
        return G1.Adj_lst == G2.Adj_lst;
    };

    //
    // operator!=
    //
    friend bool operator!=(const graph<vT, wT> &G1, const graph<vT, wT> &G2)
    {
        return G1.Adj_lst != G2.Adj_lst;
    };

    //
    // vertex_count
    //
    // Returning the # of G currently in the graph.
    //
    int vertex_count() const;

    //
    // edge_count
    //
    // Returning the # of edges currently in the graph.
    //
    int edge_count() const;

    //
    // add_vertex
    //
    // Adding the vertex v to the graph, and if so
    // returning true. If the vertex already
    // exists in the graph, then returning false.
    //
    bool add_vertex(const vT &V);

    //
    // erase_vertex
    //
    // Removing a vertex from the graph, returning true
    // if successful. If the vertex doesn't exist,
    // returning false
    //
    bool remove_vertex(const vT &V);

    //
    // add_edge
    //
    // Adding the edge (from, to, weight) to the graph, and returning
    // true.  If the Adj_lst do not exist returning false.
    //
    // NOTE: if the edge already exists, then overwriting existing
    // edge weight with the new edge weight.
    //
    // If the cycle parameter = false, the edge will not be added
    // i.e. false will be returned if the added edge forms a cycle
    //
    bool add_edge(const vT &from, const vT &to, const wT &weight = wT(), const bool &cycle = true);

    //
    // remove_edge
    //
    // Removing an edge between two Adj_lst uniquely defined by its weight
    // and returning true. Returning false if the edge does not exist
    // or at least one of the Adj_lst doesn't exist
    //
    bool remove_edge(const vT &from, const vT &to);
    
    //
    // contains_vertex
    //
    // Returning true if the function argument is present in the graph
    // else returning false
    //
    bool contains_vertex(const vT &V) const;
    
    //
    // contains_edge
    //
    // Returning true if the edge between function arguments from and to
    // exists else returning false
    //
    bool contains_edge(const vT &from, const vT &to) const;

    //
    // get_weight
    //
    // Returning the weight associated with a given edge.
    // If either of the Adj_lst don't exist or an edge between them
    // does not exist, throwing a runtime error.
    // Else returning the weight associated with that edge
    //
    wT get_weight(const vT &from, const vT &to) const;

    //
    // in_deg
    //
    // Returning the in-degree of a vertex
    //
    int in_deg(const vT &V) const;

    //
    // out_deg
    //
    // Returning the out-degree of a vertex
    //
    int out_deg(const vT &V) const;

    //
    // neighbors
    //
    // Returning a set containing the neighbors of v, i.e. all
    // Adj_lst that can be reached from v along onne edge.
    //
    set<vT> neighbors(const vT &V) const;

    //
    // get_vertices
    //
    // Returning a set containing all the Adj_lst currently in
    // the graph.
    //
    set<vT> get_vertices() const;

    //
    // connected
    //
    // checking whether a graph is connected i.e. whether
    // all vertices are reachable from every other vertex
    // used in Kruskal's algorithm
    //
    bool connected();

    //
    // BFS
    //
    // Implementing breadth-first search
    //
    // does not traverse all vertices for a disconnected graph
    //
    friend vector<vT> BFS(const graph<vT, wT> &G, const vT &start)
    {
        if (not G.Adj_lst.count(start)) {
            return {};
        }

        vector<vT> visited;
        set<vT> discovered;
        queue<vT> frontier;

        discovered.insert(start);
        frontier.push(start);

        while(not frontier.empty()) {
            vT curr = frontier.front();
            frontier.pop();

            visited.push_back(curr);

            for (const auto &adj : G.neighbors(curr)) {
                if (not discovered.count(adj)) {
                    frontier.push(adj);
                    discovered.insert(adj);
                }
            }
        }
        return visited;
    };
    //
    // overloading BFS to take a function argument and return a vector
    // of the outputs after performing any operation on the elements
    //
    // does not traverse all vertices for a disconnected graph
    //
    template<typename rT>
    friend vector<rT> BFS(const graph<vT, wT> &G, const vT &start,
                          rT(*func)(vT) /*Function pointer*/)
    {
        if (not G.Adj_lst.count(start)) {
            return {};
        }

        vector<rT> visited;
        set<vT> discovered;
        queue<vT> frontier;

        discovered.insert(start);
        frontier.push(start);

        while(not frontier.empty()) {
            vT curr = frontier.front();
            frontier.pop();

            visited.push_back(func(curr));

            for (const auto &adj : G.neighbors(curr)) {
                if (not discovered.count(adj)) {
                    frontier.push(adj);
                    discovered.insert(adj);
                }
            }
        }

        return visited;
    };

    //
    // DFS
    //
    // Implementing depth-first search
    // Also possible to implement by just replacing queue with
    // stack in the implementation of BFS
    //
    // does not traverse all vertices for a disconnected graph
    // used to check whether the graph is disconnected or not
    //
    friend vector<vT> DFS(const graph<vT, wT> &G, const vT &start)
    {
        if (not G.Adj_lst.count(start)) {
            return {};
        }

        vector<vT> visited;
        set<vT> discovered;
        stack<vT> frontier;

        frontier.push(start);

        while(not frontier.empty()) {
            vT curr = frontier.top();
            frontier.pop();
            if (not discovered.count(curr)) {
                discovered.insert(curr);

                visited.push_back(curr);

                for (const auto &adj : G.neighbors(curr)) {
                    frontier.push(adj);
                }
            }
        }

        return visited;
    };
    //
    // overloading DFS to take a function argument and return a vector
    // of the outputs after performing an operation on the elements
    //
    // does not traverse all vertices for a disconnected graph
    // used to check whether the graph is disconnected or not
    //
    template<typename rT>
    friend vector<rT> DFS(const graph<vT, wT> &G, const vT &start,
                          rT(*func)(vT) /*Function pointer*/)
    {
        if (not G.Adj_lst.count(start)) {
            return {};
        }

        vector<vT> visited;
        set<vT> discovered;
        stack<vT> frontier;

        frontier.push(start);

        while(not frontier.empty()) {
            vT curr = frontier.top();
            frontier.pop();
            if (not discovered.count(curr)) {
                discovered.insert(curr);

                visited.push_back(func(curr));

                for (const auto &adj : G.neighbors(curr)) {
                    frontier.push(adj);
                }
            }
        }

        return visited;
    };

    //
    // operator<<
    //
    // Printing the graph to an output stream
    //
    friend ostream &operator<<(ostream &out, const graph<vT, wT> &G)
    {
        out << "***************************************************" << endl;
        out << "********************* GRAPH ***********************" << endl;
        out << "**Vertex count: " << G.Adj_lst.size() << endl;
        out << "**Edge count: " << G.edge_count() << endl;
        out << "**Vertices: ";
        bool comma = false;
        for (const auto &V : G.get_vertices()) {
            if (not comma) {
                out << V;
                comma = true;
            } else {
                out << ", " << V;
            }
        }
        out << endl;
        out << "**Edges:" << endl;
        bool bar = true;
        for (const auto &V : G.get_vertices()) {
            bar = true;
            for (const auto &e : G.Adj_lst.at(V)) {
                if (bar) {
                    out << "| (" << V << ", " << e.first << ") : " << e.second << " |";
                    bar = false;
                } else
                    out << " (" << V << ", " << e.first << ") : " << e.second << " |";
            }
            out << endl;
        }
        out << "**************************************************" << endl;
        out << "**************************************************" << endl;
        return out;
    };

    //
    // MST
    //
    // Implementing Kruskal's algorithm and returning
    // the minimum spanning tree
    //
    // The concept of a spanning tree is undefined
    // for a disconnected graph
    //
    // AN MST CANNOT BE A DIGRAPH
    // Hence returning an undirected graph such that
    // when an edge from A to B is added, an edge from
    // B to A with the same weight is added
    // irrespective of whether it is present in the original graph
    // Hence onus is upon the user to find the MST of an
    // undirected graph in order not to be cheated!
    //
    graph<vT, wT> MST();
};


















//
// class multigraph
//
// graph with multiple edges between two vertices
// an edge is uniquely defined by its starting and ending
// vertices and its weight
//
template<typename vT, typename wT>
class multigraph {
private:
    //
    // Adjacency list implementation
    //
    unordered_map<vT, map<vT, set<wT>>> Adj_lst;

    //
    // forms_cycle
    //
    // Private member helper function
    // Using in add_edge to check whether an added edge
    // forms a cycle
    //
    bool forms_cycle(const vT &from, const vT &to, deque<vT> cycle, set<vT> discovered);

    //
    // connected
    //
    // Private member helper function
    // checking whether a graph is connected i.e. whether
    // all vertices are reachable from every other vertex
    // used in Kruskal's algorithm
    //
    bool connected();

public:
    //
    // default constructor
    //
    multigraph() = default;

    //
    // adding vertices for a graph from a set
    //
    multigraph(const set<vT> &vertices);

    //
    // initializing a multigraph from a graph
    //
    multigraph(const graph<vT, wT> &G);

    //
    // operator==
    //
    friend bool operator==(const multigraph<vT, wT> &G1, const multigraph<vT, wT> &G2)
    {
        return G1.Adj_lst == G2.Adj_lst;
    };

    //
    // operator!=
    //
    friend bool operator!=(const multigraph<vT, wT> &G1, const multigraph<vT, wT> &G2)
    {
        return G1.Adj_lst != G2.Adj_lst;
    };

    //
    // vertex_count
    //
    // Returning the # of G currently in the graph.
    //
    int vertex_count() const;

    //
    // edge_count
    //
    // Returning the # of edges currently in the graph.
    //
    int edge_count() const;

    //
    // add_vertex
    //
    // Adding the vertex v to the graph, and if so
    // returning true. If the vertex already
    // exists in the graph, then returning false.
    //
    bool add_vertex(const vT &V);

    //
    // erase_vertex
    //
    // Removing a vertex from the graph, returning true
    // if successful. If the vertex doesn't exist,
    // returning false
    //
    bool remove_vertex(const vT &V);

    //
    // add_edge
    //
    // Adding the edge (from, to, weight) to the graph, and returning
    // true.  If the Adj_lst do not exist returning false.
    //
    // If the cycle parameter, the edge will not be added i.e. false will
    // be returned if the added edge forms a cycle
    //
    bool add_edge(const vT &from, const vT &to, const wT &weight = wT(), const bool &cycle = true);

    //
    // remove_edge
    //
    // Removing an edge between two Adj_lst uniquely defined by its weight
    // and returning true. Returning false if the edge does not exist
    // or at least one of the Adj_lst doesn't exist
    //
    bool remove_edge(const vT &from, const vT &to, const wT &weight);

    //
    // contains_vertex
    //
    // Returning true if the function argument is present in the graph
    // else returning false
    //
    bool contains_vertex(const vT &V) const;

    //
    // contains_edge
    //
    // Returning true if the edge between function arguments from and to
    // exists else returning false
    //
    bool contains_edge(const vT &from, const vT &to, const wT &weight) const;

    //
    // num_edges
    //
    // Returning the number of edges between two vertices
    //
    int num_edges(const vT &from, const vT &to) const;

    //
    // get_weight
    //
    // Returning the weight associated with a given edge.
    // If either of the Adj_lst don't exist or an edge between them
    // does not exist, throwing a runtime error.
    // Else returning the weight associated with that edge
    //
    set<wT> get_weights(const vT &from, const vT &to) const;

    //
    // in_deg
    //
    // Returning the in-degree of a vertex
    //
    int in_deg(const vT &V) const;

    //
    // out_deg
    //
    // Returning the out-degree of a vertex
    //
    int out_deg(const vT &V) const;

    //
    // neighbors
    //
    // Returning a set containing the neighbors of v, i.e. all
    // Adj_lst that can be reached from v along onne edge.
    //
    set<vT> neighbors(const vT &V) const;

    //
    // get_vertices
    //
    // Returning a set containing all the Adj_lst currently in
    // the graph.
    //
    set<vT> get_vertices() const;

    //
    // BFS
    //
    // Implementing breadth-first search
    //
    // does not traverse all vertices for a disconnected graph
    //
    friend vector<vT> BFS(const multigraph<vT, wT> &G, const vT &start)
    {
        if (not G.Adj_lst.count(start)) {
            return {};
        }

        vector<vT> visited;
        set<vT> discovered;
        queue<vT> frontier;

        discovered.insert(start);
        frontier.push(start);

        while(not frontier.empty()) {
            vT curr = frontier.front();
            frontier.pop();

            visited.push_back(curr);

            for (const auto &adj : G.neighbors(curr)) {
                if (not discovered.count(adj)) {
                    frontier.push(adj);
                    discovered.insert(adj);
                }
            }
        }
        return visited;
    };
    //
    // overloading BFS to take a function argument and return a vector
    // of the outputs after performing any operation on the elements
    //
    // does not traverse all vertices for a disconnected graph
    //
    template<typename rT>
    friend vector<rT> BFS(const multigraph<vT, wT> &G, const vT &start,
                          rT(*func)(vT) /*Function pointer*/)
    {
        if (not G.Adj_lst.count(start)) {
            return {};
        }

        vector<rT> visited;
        set<vT> discovered;
        queue<vT> frontier;

        discovered.insert(start);
        frontier.push(start);

        while(not frontier.empty()) {
            vT curr = frontier.front();
            frontier.pop();

            visited.push_back(func(curr));

            for (const auto &adj : G.neighbors(curr)) {
                if (not discovered.count(adj)) {
                    frontier.push(adj);
                    discovered.insert(adj);
                }
            }
        }

        return visited;
    };


    //
    // DFS
    //
    // Implementing depth-first search
    // Also possible to implement by just replacing queue with
    // stack in the implementation of BFS
    //
    // does not traverse all vertices for a disconnected graph
    // used to check whether the graph is disconnected or not
    //
    friend vector<vT> DFS(const multigraph<vT, wT> &G, const vT &start)
    {
        if (not G.Adj_lst.count(start)) {
            return {};
        }

        vector<vT> visited;
        set<vT> discovered;
        stack<vT> frontier;

        frontier.push(start);

        while(not frontier.empty()) {
            vT curr = frontier.top();
            frontier.pop();
            if (not discovered.count(curr)) {
                discovered.insert(curr);

                visited.push_back(curr);

                for (const auto &adj : G.neighbors(curr)) {
                    frontier.push(adj);
                }
            }
        }

        return visited;
    };
    //
    // overloading DFS to take a function argument and return a vector
    // of the outputs after performing an operation on the elements
    //
    // does not traverse all vertices for a disconnected graph
    //
    template<typename rT>
    friend vector<rT> DFS(const multigraph<vT, wT> &G, const vT &start,
                          rT(*func)(vT) /*Function pointer*/)
    {
        if (not G.Adj_lst.count(start)) {
            return {};
        }

        vector<vT> visited;
        set<vT> discovered;
        stack<vT> frontier;

        frontier.push(start);

        while(not frontier.empty()) {
            vT curr = frontier.top();
            frontier.pop();
            if (not discovered.count(curr)) {
                discovered.insert(curr);

                visited.push_back(func(curr));

                for (const auto &adj : G.neighbors(curr)) {
                    frontier.push(adj);
                }
            }
        }

        return visited;
    };

    //
    // operator<<
    //
    // Printing the graph to an output stream
    //
    friend ostream &operator<<(ostream &out, const multigraph<vT, wT> &G)
    {
        out << "***************************************************" << endl;
        out << "********************* GRAPH ***********************" << endl;
        out << "**Vertex count: " << G.Adj_lst.size() << endl;
        out << "**Edge count: " << G.edge_count() << endl;
        out << "**Vertices: ";
        bool comma = false;
        for (const auto &V : G.get_vertices()) {
            if (not comma) {
                out << V;
                comma = true;
            } else {
                out << ", " << V;
            }
        }
        out << endl;
        out << "**Edges:" << endl;
        bool bar = true;
        for (const auto &V : G.get_vertices()) {
            bar = true;
            for (const auto &e : G.Adj_lst.at(V)) {
                for (const auto &W : e.second) {
                    if (bar) {
                        out << "| (" << V << ", " << e.first << ") : " << W << " |";
                        bar = false;
                    } else
                        out << " (" << V << ", " << e.first << ") : " << W << " |";
                }
            }
            out << endl;
        }
        out << "**************************************************" << endl;
        out << "**************************************************" << endl;
        return out;
    };

    //
    // MST
    //
    // Implementing Kruskal's algorithm and returning
    // the minimum spanning tree
    //
    // The concept of a spanning tree is undefined
    // for a disconnected graph
    //
    // AN MST CANNOT BE A DIGRAPH
    // Hence returning an undirected graph such that
    // when an edge from A to B is added, an edge from
    // B to A with the same weight is added
    // irrespective of whether it is present in the original graph
    // Hence onus is upon the user to find the MST of an
    // undirected graph in order not to be cheated!
    //
    graph<vT, wT> MST();
};












































































































//
// vv IMPLEMENTATIONS ABSTRACTED AWAY vv
//

































































































//
// FUNCTION DEFINITIONS
//
template<typename vT, typename wT>
bool graph<vT, wT>::forms_cycle(const vT &from, const vT &to, deque<vT> cycle, set<vT> discovered)
{
    discovered.insert(from);
    cycle.push_back(from);
    for (const auto &adj : neighbors(to)) {
        if (adj == cycle.front() or (not discovered.count(adj) and
            forms_cycle(to, adj, cycle, discovered))) {
            return true;
        }
    }
    return false;
}

template<typename vT, typename wT>
bool graph<vT, wT>::connected()
{
    set<vT> traversed;
    set<vT> vertices = get_vertices();
    for (const auto &V : DFS(*this, *vertices.begin())) {
        traversed.insert(V);
    }
    return traversed == vertices;
}

template<typename vT, typename wT>
int graph<vT, wT>::vertex_count() const
{
    return this->Adj_lst.size();
}

template<typename vT, typename wT>
int graph<vT, wT>::edge_count() const
{
    int ecount = 0;
    for (const auto &e : Adj_lst) {
        ecount += e.second.size();
    }
    return ecount;
}

template<typename vT, typename wT>
bool graph<vT, wT>::add_vertex(const vT &V)
{
    // is the vertex already in the graph?  If so, we do not
    // insert again otherwise Vertices may fill with duplicates:
    if (Adj_lst.count(V)) {
        return false;
    }
    // if here is reached, vertex does not exist so insert.
    Adj_lst[V];
    return true;
}

template<typename vT, typename wT>
graph<vT, wT>::graph(const set<vT> &vertices)
{
    for (const auto &V : vertices) {
        add_vertex(V);
    }
}

template<typename vT, typename wT>
bool graph<vT, wT>::remove_vertex(const vT &V)
{
    // is the vertex not in the graph?
    if (not Adj_lst.count(V)) {
        return false;
    }
    // if here is reached, then the vertex is present
    for (auto &e : Adj_lst) {
        if (e.second.count(V)) e.second.erase(V);
    }
    Adj_lst.erase(V);
    return true;
}

template<typename vT, typename wT>
bool graph<vT, wT>::add_edge(const vT &from, const vT &to, const wT &weight, const bool &cycle)
{
    // Are both the Adj_lst in the graph?
    if (not Adj_lst.count(from) or not Adj_lst.count(to)) {
        return false;
    }
    // If here is reached then they are present
    Adj_lst.at(from)[to] = weight;
    if (cycle or not forms_cycle(from, to, {}, {})) {
        return true;
    } else  {
        Adj_lst.at(from).erase(to);
        return false;
    }
}

template<typename vT, typename wT>
bool graph<vT, wT>::remove_edge(const vT &from, const vT &to)
{
    // Are both Adj_lst in the graph?
    // Does an edge exist between them if they are present?
    if (not Adj_lst.count(from) or not Adj_lst.count(to) or not Adj_lst.at(from).count(to)) {
        return false;
    }
    // If here is reached, an edge does exist between them
    Adj_lst.at(from).erase(to);
    return true;
}

template<typename vT, typename wT>
bool graph<vT, wT>::contains_vertex(const vT &V) const
{
    return Adj_lst.count(V);
}

template<typename vT, typename wT>
bool graph<vT, wT>::contains_edge(const vT &from, const vT &to) const
{
    return Adj_lst.count(from) and Adj_lst.count(to) and Adj_lst.at(from).count(to);
}

template<typename vT, typename wT>
wT graph<vT, wT>::get_weight(const vT &from, const vT &to) const
{
    if (not contains_edge(from, to)) {
        throw runtime_error("/*** get_weight: AT LEAST ONE VERTEX OR EDGE DOES NOT EXIST ***/");
    }
    return Adj_lst.at(from).at(to);
}

template<typename vT, typename wT>
int graph<vT, wT>::in_deg(const vT &V) const
{
    if (not Adj_lst.count(V)) {
        return -1;
    }
    int in = 0;
    for (const auto &e : Adj_lst) {
        if (e.second.count(V)) ++in;
    }
    return in;
}

template<typename vT, typename wT>
int graph<vT, wT>::out_deg(const vT &V) const
{
    if (not Adj_lst.count(V)) {
        return -1;
    }
    return Adj_lst.at(V).size();
}

template<typename vT, typename wT>
set<vT> graph<vT, wT>::neighbors(const vT &V) const
{
    set<vT> neighbors;
    if (Adj_lst.count(V)) {
        for (const auto &e : Adj_lst.at(V)) {
            neighbors.insert(e.first);
        }
    }
    return neighbors;
}

template<typename vT, typename wT>
set<vT> graph<vT, wT>::get_vertices() const
{
    set<vT> vertices;
    for (const auto &e : Adj_lst) {
        vertices.insert(e.first);
    }
    return vertices;
}

template<typename vT, typename wT>
graph<vT, wT> graph<vT, wT>::MST()
{
    if (not connected()) {
        throw runtime_error("/*** MST: UNDEFINED FOR A DISCONNECTED GRAPH ***/");
    }

    graph<vT, wT> mst;
    class prioritize {
    public:
        bool operator()(const pair<pair<vT, vT>, wT> &p1, const pair<pair<vT, vT>, wT> &p2)
        {
            return p1.second > p2.second;
        }
    };
    priority_queue<pair<pair<vT, vT>, wT>,
            vector<pair<pair<vT, vT>, wT>>,
            prioritize> edges;

    for (const auto &V : get_vertices()) {
        mst.add_vertex(V);
        for (const auto &e : Adj_lst.at(V)) {
            edges.push({{V, e.first}, e.second});
        }
    }

    while (not mst.connected()) {
        if (not mst.contains_edge(edges.top().first.second, edges.top().first.first) and
            not mst.contains_edge(edges.top().first.first, edges.top().first.second)) {
            // Adding an edge from A to B
            if (mst.add_edge(edges.top().first.first, edges.top().first.second,
                             edges.top().second, false)) {

                // Adding an edge from B to A with same weight
                // Hence ensuring undirectedness of the MST
                mst.add_edge(edges.top().first.second, edges.top().first.first,
                             edges.top().second);
            }
        }
        edges.pop();
    }

    return mst;
}


//
// FUNCTION DEFINITIONS for MULTIGRAPH
//
template<typename vT, typename wT>
bool multigraph<vT, wT>::forms_cycle(const vT &from, const vT &to, deque<vT> cycle, set<vT> discovered)
{
    discovered.insert(from);
    cycle.push_back(from);
    for (const auto &adj : neighbors(to)) {
        if (adj == cycle.front() or (not discovered.count(adj) and
                                     forms_cycle(to, adj, cycle, discovered))) {
            return true;
        }
    }
    return false;
}

template<typename vT, typename wT>
bool multigraph<vT, wT>::connected()
{
    set<vT> traversed;
    set<vT> vertices = get_vertices();
    for (const auto &V : DFS(*this, *vertices.begin())) {
        traversed.insert(V);
    }
    return traversed == vertices;
}

template<typename vT, typename wT>
int multigraph<vT, wT>::vertex_count() const
{
    return this->Adj_lst.size();
}

template<typename vT, typename wT>
int multigraph<vT, wT>::edge_count() const
{
    int ecount = 0;
    for (const auto &V : Adj_lst) {
        for (const auto &e : V.second) {
            ecount += e.second.size();
        }
    }
    return ecount;
}

template<typename vT, typename wT>
bool multigraph<vT, wT>::add_vertex(const vT &V)
{
    // is the vertex already in the graph?  If so, we do not
    // insert again otherwise Vertices may fill with duplicates:
    if (Adj_lst.count(V)) {
        return false;
    }
    // if here is reached, vertex does not exist so insert.
    Adj_lst[V];
    return true;
}

template<typename vT, typename wT>
multigraph<vT, wT>::multigraph(const set<vT> &vertices)
{
    for (const auto &V : vertices) {
        add_vertex(V);
    }
}

template<typename vT, typename wT>
multigraph<vT, wT>::multigraph(const graph<vT, wT> &G)
{
    set<vT> vertices = G.get_vertices();
    for (const auto &from : vertices) {
        add_vertex(from);
        for (const auto &to : vertices) {
            if (G.contains_edge(from, to)) {
                add_edge(from, to, G.get_weight(from, to));
            }
        }
    }
}

template<typename vT, typename wT>
bool multigraph<vT, wT>::remove_vertex(const vT &V)
{
    // is the vertex not in the graph?
    if (not Adj_lst.count(V)) {
        return false;
    }
    // if here is reached, then the vertex is present
    for (auto &e : Adj_lst) {
        if (e.second.count(V)) e.second.erase(V);
    }
    Adj_lst.erase(V);
    return true;
}

template<typename vT, typename wT>
bool multigraph<vT, wT>::add_edge(const vT &from, const vT &to, const wT &weight, const bool &cycle)
{
    // Are both the Adj_lst in the graph?
    if (not Adj_lst.count(from) or not Adj_lst.count(to)) {
        return false;
    }
    // If here is reached then they are present
    Adj_lst.at(from)[to].insert(weight);
    if (cycle or not forms_cycle(from, to, {}, {})) {
        return true;
    } else  {
        Adj_lst.at(from).at(to).erase(weight);
        return false;
    }
}

template<typename vT, typename wT>
bool multigraph<vT, wT>::remove_edge(const vT &from, const vT &to, const wT &weight)
{
    // Are both Adj_lst in the graph?
    // Does an edge exist between them if they are present?
    if (not Adj_lst.count(from) or not Adj_lst.count(to) or not Adj_lst.at(from).count(to)) {
        return false;
    }
    // If here is reached, an edge does exist between them
    Adj_lst.at(from).at(to).erase(weight);
    return true;
}

template<typename vT, typename wT>
bool multigraph<vT, wT>::contains_vertex(const vT &V) const
{
    return Adj_lst.count(V);
}

template<typename vT, typename wT>
bool multigraph<vT, wT>::contains_edge(const vT &from, const vT &to, const wT &weight) const
{
    return Adj_lst.count(from) and Adj_lst.count(to) and Adj_lst.at(from).at(to).count(weight);
}

template<typename vT, typename wT>
int multigraph<vT, wT>::num_edges(const vT &from, const vT &to) const
{
    return Adj_lst.at(from).at(to).size();
}

template<typename vT, typename wT>
set<wT> multigraph<vT, wT>::get_weights(const vT &from, const vT &to) const
{
    if (not Adj_lst.count(from) or not Adj_lst.count(to) or not num_edges(from, to)) {
        throw runtime_error("/*** get_weight: AT LEAST ONE VERTEX OR EDGES DO NOT EXIST ***/");
    }
    return Adj_lst.at(from).at(to);
}

template<typename vT, typename wT>
int multigraph<vT, wT>::in_deg(const vT &V) const
{
    if (not Adj_lst.count(V)) {
        return -1;
    }
    int in = 0;
    for (const auto &e : Adj_lst) {
        if (e.second.count(V))  in += e.second.second.size();
    }
    return in;
}

template<typename vT, typename wT>
int multigraph<vT, wT>::out_deg(const vT &V) const
{
    if (not Adj_lst.count(V)) {
        return -1;
    }
    int out = 0;
    for (const auto &e : Adj_lst.at(V)) {
        out += e.second.size();
    }
    return out;
}

template<typename vT, typename wT>
set<vT> multigraph<vT, wT>::neighbors(const vT &V) const
{
    set<vT> neighbors;
    if (Adj_lst.count(V)) {
        for (const auto &e : Adj_lst.at(V)) {
            neighbors.insert(e.first);
        }
    }
    return neighbors;
}

template<typename vT, typename wT>
set<vT> multigraph<vT, wT>::get_vertices() const
{
    set<vT> vertices;
    for (const auto &e : Adj_lst) {
        vertices.insert(e.first);
    }
    return vertices;
}

template<typename vT, typename wT>
graph<vT, wT> multigraph<vT, wT>::MST()
{
    if (not connected()) {
        throw runtime_error("/*** MST: UNDEFINED FOR A DISCONNECTED GRAPH ***/");
    }

    graph<vT, wT> mst;
    class prioritize {
    public:
        bool operator()(const pair<pair<vT, vT>, wT> &p1, const pair<pair<vT, vT>, wT> &p2)
        {
            return p1.second > p2.second;
        }
    };
    priority_queue<pair<pair<vT, vT>, wT>,
            vector<pair<pair<vT, vT>, wT>>,
            prioritize> edges;

    for (const auto &V : get_vertices()) {
        mst.add_vertex(V);
        for (const auto &e : Adj_lst.at(V)) {
            for (const auto &w : e.second) {
                edges.push({{V, e.first}, w});
            }
        }
    }

    while (not mst.connected()) {
        if (not mst.contains_edge(edges.top().first.second, edges.top().first.first) and
            not mst.contains_edge(edges.top().first.first, edges.top().first.second)) {
            // Adding an edge from A to B
            if (mst.add_edge(edges.top().first.first, edges.top().first.second,
                             edges.top().second, false)) {

                // Adding an edge from B to A with same weight
                // Hence ensuring undirectedness of the MST
                mst.add_edge(edges.top().first.second, edges.top().first.first,
                             edges.top().second);
            }
        }
        edges.pop();
    }

    return mst;
}