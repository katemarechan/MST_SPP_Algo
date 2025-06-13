#ifndef MST_SPP_HPP
#define MST_SPP_HPP

#include <iostream>
#include <vector>
#include <queue>
#include <climits>
#include <ctime>
#include <chrono>
#include <random>
#include <algorithm>
#include <tuple>

namespace MST_SPP {

    using namespace std;

    // Structure to represent an edge with a weight
    struct Edge {
        int src, dest, weight;
    };

    // Adjacency Matrix class to represent the graph
    class AdjacencyMatrix {
    public:
        // Constructor to initialize an adjacency matrix for V vertices
        AdjacencyMatrix(int V) : V(V) {
            adj = new int* [V];
            for (int i = 0; i < V; i++) {
                adj[i] = new int[V];
                for (int j = 0; j < V; j++) {
                    adj[i][j] = 0;
                }
            }
        }

        // Destructor to deallocate the adjacency matrix
        ~AdjacencyMatrix() {
            for (int i = 0; i < V; i++) {
                delete[] adj[i];
            }
            delete[] adj;
        }

        // Function to add an edge to the adjacency matrix
        void addEdge(int src, int dest, int weight) {
            adj[src][dest] = weight;
            adj[dest][src] = weight;
        }

        // Function to get the adjacency matrix (const version)
        int** getMatrix() const {
            return adj;
        }

        // Function to get the adjacency matrix (non-const version)
        int** getMatrix() {
            return adj;
        }

        // Function to get number of vertices (const version)
        int getV() const {
            return V;
        }

        // Function to get number of vertices (non-const version)
        int& getV() {
            return V;
        }

    private:
        int V;
        int** adj;
    };

    class AdjacencyList {
    public:
        // Constructor to initialize an adjacency list for V vertices
        AdjacencyList(int V) : V(V) {
            adj = new vector<Edge>[V];
        }

        // Destructor to deallocate the adjacency list
        ~AdjacencyList() {
            delete[] adj;
        }

        // Function to add an edge to the adjacency list
        void addEdge(int src, int dest, int weight) {
            // Helper lambda to check if edge already exists
            auto edgeExists = [](const vector<Edge>& edgeList, int destination) {
                for (const auto& edge : edgeList) {
                    if (edge.dest == destination) {
                        return true;
                    }
                }
                return false;
            };

            // Only add edge from src to dest if it doesn't already exist
            if (!edgeExists(adj[src], dest)) {
                adj[src].push_back({ src, dest, weight });
            }

            // Only add edge from dest to src if it doesn't already exist
            if (!edgeExists(adj[dest], src)) {
                adj[dest].push_back({ dest, src, weight });
            }
        }

        // Function to get the adjacency list (const version)
        const vector<Edge>* getList() const {
            return adj;
        }

        // Function to get the adjacency list (non-const version)
        vector<Edge>* getList() {
            return adj;
        }

        // Function to get number of vertices (const version)
        int getV() const {
            return V;
        }

        // Function to get number of vertices (non-const version)
        int& getV() {
            return V;
        }

    private:
        int V;
        vector<Edge>* adj;
    };

    // Function declarations - implementations will be in separate .cpp file
    std::vector<Edge> primMSTMatrix(AdjacencyMatrix& graph);
    std::vector<Edge> primMSTList(AdjacencyList& graph);
    int findRoot(vector<int>& parent, int i);
    void unionByRank(vector<int>& parent, vector<int>& rank, int x, int y);
    std::vector<Edge> kruskalMSTMatrix(AdjacencyMatrix& graph);
    vector<Edge> kruskalMSTList(vector<Edge>* adj, int V);
    std::pair<std::vector<int>, std::vector<int>> dijkstraList(AdjacencyList& graph, int start);
    std::pair<std::vector<int>, std::vector<int>> dijkstraMatrix(int** adj, int V, int start);
    std::pair<std::vector<int>, std::vector<int>> FordBellmanMatrix(AdjacencyMatrix& graph, int start);
    std::pair<std::vector<int>, std::vector<int>> FordBellmanList(std::vector<Edge>* adj, int V, int start);
    int calculateTotalWeight(const vector<Edge>& mst);
    void printMSTforMatrix(const std::vector<Edge>& mst, const AdjacencyMatrix& graph);
    void printMSTforList(const vector<Edge>& mst);
    void displayShortestPathsWithPredecessors(const std::vector<int>& distances, const std::vector<int>& predecessors, int source);
    void testAdjencyMatrixMST(AdjacencyMatrix& obj, int n);
    void testAdjencyMatrixSPP(AdjacencyMatrix& obj, int n, int dijkstraStart);
    void testAdjencyListMST(AdjacencyList& obj, int n);
    void testAdjencyListSPP(AdjacencyList& obj, int n, int dijkstraStart);

} // namespace MST_SPP

#endif // MST_SPP_HPP