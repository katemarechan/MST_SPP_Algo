#include "MST_SPP.hpp"

namespace MST_SPP {

    /**
     * PRIM'S ALGORITHM - MATRIX IMPLEMENTATION
     * Algorithm:
     * 1. Start with arbitrary vertex (vertex 0)
     * 2. Mark it as visited and add all its edges to priority queue
     * 3. While queue not empty:
     *    - Extract minimum weight edge
     *    - If destination not visited, add edge to MST
     *    - Mark destination as visited
     *    - Add all edges from new vertex to queue
     */
    std::vector<Edge> primMSTMatrix(AdjacencyMatrix& graph) {
        std::vector<Edge> mst;                                        // Store MST edges
        std::vector<bool> visited(graph.getV(), false);              // Track visited vertices
        // Priority queue stores: (negative_weight, destination_vertex, source_vertex)
        std::priority_queue<std::tuple<int, int, int>> pq;

        // INITIALIZATION PHASE: Start from vertex 0
        visited[0] = true;                                            // Mark starting vertex as visited

        // Add all edges from vertex 0 to priority queue
        for (int i = 0; i < graph.getV(); i++) {
            if (graph.getMatrix()[0][i] != 0) {                       // If edge exists
                // Store: (negative_weight, destination, source)
                pq.push(std::make_tuple(-graph.getMatrix()[0][i], i, 0));
            }
        }

        // MAIN ALGORITHM LOOP
        while (!pq.empty() && mst.size() < graph.getV() - 1) {
            int weight = -std::get<0>(pq.top());                      // Get actual weight (un-negate)
            int dest = std::get<1>(pq.top());                         // Get destination vertex
            int src = std::get<2>(pq.top());                          // Get source vertex
            pq.pop();                                                 // Remove from queue

            if (!visited[dest]) {                                     // If destination not yet in MST
                visited[dest] = true;                                 // Mark as visited
                mst.push_back({src, dest, weight});                   // Add edge to MST

                // Add all edges from newly added vertex to queue
                for (int v = 0; v < graph.getV(); v++) {
                    if (graph.getMatrix()[dest][v] != 0 && !visited[v]) { // If edge exists and destination not visited
                        pq.push(std::make_tuple(-graph.getMatrix()[dest][v], v, dest));
                    }
                }
            }
        }

        return mst;                                                   // Return minimum spanning tree
    }

    /**
     * PRIM'S ALGORITHM - ADJACENCY LIST IMPLEMENTATION
     * More efficient for sparse graphs due to better edge iteration
     */
    std::vector<Edge> primMSTList(AdjacencyList& graph) {
        std::vector<Edge> mst;                                        // Store MST edges
        std::vector<bool> visited(graph.getV(), false);              // Track visited vertices
        // Priority queue stores: (negative_weight, destination_vertex, source_vertex)
        std::priority_queue<std::tuple<int, int, int>> pq;

        // INITIALIZATION PHASE: Start from vertex 0
        visited[0] = true;                                            // Mark starting vertex as visited

        // Add all edges from vertex 0 to priority queue
        for (const auto& e : graph.getList()[0]) {                   // Iterate through adjacency list
            // Store: (negative_weight, destination, source)
            pq.push(std::make_tuple(-e.weight, e.dest, 0));
        }

        // MAIN ALGORITHM LOOP
        while (!pq.empty() && mst.size() < graph.getV() - 1) {
            int weight = -std::get<0>(pq.top());                      // Get actual weight (un-negate)
            int dest = std::get<1>(pq.top());                         // Get destination vertex
            int src = std::get<2>(pq.top());                          // Get source vertex
            pq.pop();                                                 // Remove from queue

            if (!visited[dest]) {                                     // If destination not yet in MST
                visited[dest] = true;                                 // Mark as visited
                mst.push_back({src, dest, weight});                   // Add edge to MST

                // Add all edges from newly added vertex to queue
                for (const auto& e : graph.getList()[dest]) {         // Iterate through adjacency list
                    if (!visited[e.dest]) {                           // If destination not visited
                        pq.push(std::make_tuple(-e.weight, e.dest, dest));
                    }
                }
            }
        }

        return mst;                                                   // Return minimum spanning tree
    }

    /**
     * UNION-FIND DATA STRUCTURE HELPER FUNCTION
     * Used in Kruskal's algorithm for cycle detection
     *
     * Path compression optimization: Makes faster for future operations
     */
    int findRoot(vector<int>& parent, int i) {
        if (parent[i] != i) {                                     // If not root of its set
            parent[i] = findRoot(parent, parent[i]);        // Path compression: point directly to root
        }
        return parent[i];                                        // Return root of the set
    }

    /**
     * UNION BY RANK OPERATION
     * Merges two sets by rank to keep the tree balanced
     * Always attaches the tree with smaller rank under the tree with larger rank
     *
     * This optimization prevents trees from becoming unbalanced chains
     * Time Complexity: Nearly O(1) amortized
     */
    void unionByRank(vector<int>& parent, vector<int>& rank, int x, int y) {
        int rootX = findRoot(parent, x);                          // Find root of x's set
        int rootY = findRoot(parent, y);                          // Find root of y's set

        if (rootX != rootY) {                                     // If vertices in different sets
            // Union by rank: attach smaller rank tree under root of higher rank tree
            if (rank[rootX] < rank[rootY]) {
                parent[rootX] = rootY;                            // Attach X's tree under Y's root
            } else if (rank[rootX] > rank[rootY]) {
                parent[rootY] = rootX;                            // Attach Y's tree under X's root
            } else {
                // If ranks are equal, make one root and increment its rank
                parent[rootY] = rootX;                            // Make X the root
                rank[rootX]++;                                    // Increment X's rank
            }
        }
    }

    /**
     * KRUSKAL'S ALGORITHM - MATRIX IMPLEMENTATION
     * Algorithm:
     * 1. Extract all edges from adjacency matrix
     * 2. Sort edges by weight (ascending)
     * 3. Use Union-Find to detect cycles
     * 4. Add edge to MST if it doesn't create cycle
     */
    std::vector<Edge> kruskalMSTMatrix(AdjacencyMatrix& graph) {
        std::vector<Edge> mst;                                    // Store MST edges
        std::vector<Edge> edges;                                  // Store all graph edges
        std::vector<int> parent(graph.getV());                 // Union-Find parent array
        std::vector<int> rank(graph.getV(), 0);                  // Union-Find rank array

        // INITIALIZATION: Each vertex is its own parent (disjoint sets)
        for (int i = 0; i < graph.getV(); i++) {
            parent[i] = i;                                        // Each vertex is root of its own set
            rank[i] = 0;                                          // Initial rank is 0
        }

        // EDGE EXTRACTION: Get all edges from adjacency matrix
        for (int i = 0; i < graph.getV(); i++) {
            for (int j = i + 1; j < graph.getV(); j++) {          // Only upper triangle (undirected graph)
                if (graph.getMatrix()[i][j] != 0) {               // If edge exists
                    edges.push_back({ i, j, graph.getMatrix()[i][j] }); // Add edge to list
                }
            }
        }

        // SORTING: Sort edges by weight (ascending order)
        std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
            return a.weight < b.weight;                           // Compare by weight
        });

        // MAIN ALGORITHM: Process edges in order of increasing weight
        for (const auto& e : edges) {
            int root1 = findRoot(parent, e.src);                  // Find root of source vertex set
            int root2 = findRoot(parent, e.dest);                 // Find root of destination vertex set

            if (root1 != root2) {                                 // If vertices in different sets (no cycle)
                mst.push_back(e);                                 // Add edge to MST
                unionByRank(parent, rank, e.src, e.dest);         // Union by rank: merge sets efficiently
            }
            // If root1 == root2, adding edge would create cycle, so skip
        }

        return mst;                                               // Return minimum spanning tree
    }

    /**
     * KRUSKAL'S ALGORITHM - ADJACENCY LIST IMPLEMENTATION
     * More efficient edge extraction for sparse graphs
     */
    vector<Edge> kruskalMSTList(vector<Edge>* adj, int V) {
        vector<Edge> mst;                                         // Store MST edges
        vector<Edge> edges;                                       // Store all graph edges
        vector<int> parent(V);                                    // Union-Find parent array
        vector<int> rank(V, 0);                                   // Union-Find rank array

        // INITIALIZATION: Each vertex is its own parent
        for (int i = 0; i < V; i++) {
            parent[i] = i;                                        // Each vertex is root of its own set
            rank[i] = 0;                                          // Initial rank is 0
        }

        // EDGE EXTRACTION: Get all edges from adjacency list
        for (int i = 0; i < V; i++) {
            for (const auto& e : adj[i]) {                        // For each edge in vertex i's list
                // Only add edge once (when src < dest) to avoid duplicates
                // Since graph is undirected, each edge appears in both adjacency lists
                if (e.src < e.dest) {
                    edges.push_back(e);                           // Add edge to collection
                }
            }
        }

        // SORTING: Sort edges by weight (ascending order)
        sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
            return a.weight < b.weight;                           // Compare by weight
        });

        // MAIN ALGORITHM: Process edges in order of increasing weight
        for (const auto& e : edges) {
            int root1 = findRoot(parent, e.src);                  // Find root of source vertex set
            int root2 = findRoot(parent, e.dest);                 // Find root of destination vertex set

            if (root1 != root2) {                                 // If vertices in different sets (no cycle)
                mst.push_back(e);                                 // Add edge to MST
                unionByRank(parent, rank, e.src, e.dest);         // Union by rank: merge sets efficiently
            }
        }

        return mst;                                               // Return minimum spanning tree
    }

    /**
     * DIJKSTRA'S ALGORITHM - ADJACENCY LIST IMPLEMENTATION
     *
     * Finds shortest paths from source to all other vertices
     * Works only with non-negative edge weights
     *
     * Algorithm:
     * 1. Initialize distances to infinity, source to 0
     * 2. Use priority queue to process vertices in order of distance
     * 3. For each vertex, relax all adjacent edges
     * 4. Mark vertex as processed to avoid reprocessing
     */
    std::pair<std::vector<int>, std::vector<int>> dijkstraList(AdjacencyList& graph, int start) {
        std::vector<int> dist(graph.getV(), INT_MAX);             // Distance array: infinity initially
        std::vector<int> predecessor(graph.getV(), -1);           // ADD: Predecessor array
        std::vector<bool> visited(graph.getV(), false);           // Track processed vertices

        // INITIALIZATION
        dist[start] = 0;                                          // Distance to source is 0
        std::priority_queue<std::pair<int, int>> pq;              // Min-heap: (negative_distance, vertex)
        pq.push({ 0, start });                                    // Add source to queue

        // MAIN ALGORITHM LOOP
        while (!pq.empty()) {
            int u = pq.top().second;                              // Get vertex with minimum distance
            pq.pop();                                             // Remove from queue

            if (visited[u]) {                                     // Skip if already processed
                continue;                                         // (priority queue may have duplicates)
            }

            visited[u] = true;                                    // Mark as processed

            // EDGE RELAXATION: Update distances to adjacent vertices
            for (const auto& e : graph.getList()[u]) {            // For each neighbor
                if (!visited[e.dest] && dist[e.dest] > dist[u] + e.weight) {          // If shorter path found
                    dist[e.dest] = dist[u] + e.weight;            // Update distance
                    predecessor[e.dest] = u;                      // ADD: Track predecessor
                    pq.push({ -dist[e.dest], e.dest });          // Add to queue (negative for min-heap)
                }
            }
        }

        return std::make_pair(dist, predecessor);                 // CHANGE: Return both
    }

    /**
     * DIJKSTRA'S ALGORITHM - MATRIX IMPLEMENTATION
     * Time Complexity: O(V²) - checks all vertices for each vertex
     * Space Complexity: O(V)
     *
     */
    std::pair<std::vector<int>, std::vector<int>> dijkstraMatrix(int** adj, int V, int start) {
        std::vector<int> dist(V, INT_MAX);                        // Distance array: infinity initially
        std::vector<int> predecessor(V, -1);                      // ADD: Predecessor array
        std::vector<bool> visited(V, false);                      // Track processed vertices

        // INITIALIZATION
        dist[start] = 0;                                          // Distance to source is 0
        std::priority_queue<std::pair<int, int>> pq;              // Min-heap: (negative_distance, vertex)
        pq.push({ 0, start });                                    // Add source to queue

        // MAIN ALGORITHM LOOP
        while (!pq.empty()) {
            int u = pq.top().second;                              // Get vertex with minimum distance
            pq.pop();                                             // Remove from queue

            if (visited[u]) {                                     // Skip if already processed
                continue;
            }

            visited[u] = true;                                    // Mark as processed

            // EDGE RELAXATION: Check all possible destinations
            for (int i = 0; i < V; i++) {
                if (adj[u][i] != 0 && dist[i] > dist[u] + adj[u][i]) { // If edge exists and shorter path found
                    dist[i] = dist[u] + adj[u][i];                // Update distance
                    predecessor[i] = u;                           // ADD: Track predecessor
                    pq.push({ -dist[i], i });                     // Add to queue (negative for min-heap)
                }
            }
        }

        return std::make_pair(dist, predecessor);                 // CHANGE: Return both
    }
    /**
     * BELLMAN-FORD ALGORITHM - MATRIX IMPLEMENTATION WITH NEGATIVE CYCLE DETECTION
     * Time Complexity: O(V³) - V iterations × V² edge checks
     * Space Complexity: O(V)
     *
     * Can handle negative edge weights (unlike Dijkstra)
     * Can detect negative cycles
     *
     * Algorithm:
     * 1. Initialize distances to infinity, source to 0
     * 2. Relax all edges V-1 times (V-1 iterations)
     * 3. Each iteration guarantees shortest paths of length ? iteration number
     * 4. After V-1 iterations, all shortest paths are found
     * 5. Run one more iteration - if any distance updates, negative cycle exists
     */
    std::pair<std::vector<int>, std::vector<int>> FordBellmanMatrix(AdjacencyMatrix& graph, int start) {
        std::vector<int> dist(graph.getV(), INT_MAX);             // Distance array: infinity initially
        std::vector<int> predecessor(graph.getV(), -1);           // ADD: Predecessor array
        dist[start] = 0;                                          // Distance to source is 0

        // MAIN ALGORITHM: Relax all edges V-1 times
        for (int i = 0; i < graph.getV() - 1; i++) {             // V-1 iterations
            // Check all possible edges (u,v)
            for (int u = 0; u < graph.getV(); u++) {             // For each source vertex u
                for (int v = 0; v < graph.getV(); v++) {         // For each destination vertex v
                    // EDGE RELAXATION: If edge exists and shorter path found
                    if (graph.getMatrix()[u][v] != 0 &&           // Edge exists
                        dist[u] != INT_MAX &&                     // Source is reachable
                        dist[v] > dist[u] + graph.getMatrix()[u][v]) { // Shorter path found
                        dist[v] = dist[u] + graph.getMatrix()[u][v];   // Update distance
                        predecessor[v] = u;                       // ADD: Track predecessor
                    }
                }
            }
        }

        // NEGATIVE CYCLE DETECTION: Run one more iteration
        for (int u = 0; u < graph.getV(); u++) {                 // For each source vertex u
            for (int v = 0; v < graph.getV(); v++) {             // For each destination vertex v
                // If we can still relax an edge, negative cycle exists
                if (graph.getMatrix()[u][v] != 0 &&               // Edge exists
                    dist[u] != INT_MAX &&                         // Source is reachable
                    dist[v] > dist[u] + graph.getMatrix()[u][v]) { // Can still relax

                    std::cout << "\n*** NEGATIVE CYCLE DETECTED ***" << std::endl;
                    std::cout << "Graph contains a negative weight cycle!" << std::endl;
                    std::cout << "Shortest path distances are not well-defined." << std::endl;
                    std::cout << "Negative cycle involves edge: (" << u << " -> " << v
                              << ") with weight " << graph.getMatrix()[u][v] << std::endl;

                    // Return empty vectors to indicate negative cycle
                    return std::make_pair(std::vector<int>(), std::vector<int>());
                }
            }
        }

        return std::make_pair(dist, predecessor);                 // CHANGE: Return both
    }

    /**
     * BELLMAN-FORD ALGORITHM - ADJACENCY LIST IMPLEMENTATION WITH NEGATIVE CYCLE DETECTION
     * Time Complexity: O(VE) where V is vertices, E is edges
     * Space Complexity: O(V)
     *
     * More efficient for sparse graphs
     */
    std::pair<std::vector<int>, std::vector<int>> FordBellmanList(std::vector<Edge>* adj, int V, int start) {
        std::vector<int> dist(V, INT_MAX);                        // Distance array: infinity initially
        std::vector<int> predecessor(V, -1);                      // ADD: Predecessor array
        dist[start] = 0;                                          // Distance to source is 0

        // MAIN ALGORITHM: Relax all edges V-1 times
        for (int i = 0; i < V - 1; i++) {                        // V-1 iterations
            // Check all edges in adjacency list
            for (int j = 0; j < V; j++) {                         // For each vertex j
                for (const auto& e : adj[j]) {                    // For each edge from vertex j
                    // EDGE RELAXATION: If shorter path found
                    if (dist[e.src] != INT_MAX &&                 // Source is reachable
                        dist[e.dest] > dist[e.src] + e.weight) {  // Shorter path found
                        dist[e.dest] = dist[e.src] + e.weight;    // Update distance
                        predecessor[e.dest] = e.src;              // ADD: Track predecessor
                    }
                }
            }
        }

        // NEGATIVE CYCLE DETECTION: Run one more iteration
        for (int j = 0; j < V; j++) {                             // For each vertex j
            for (const auto& e : adj[j]) {                        // For each edge from vertex j
                // If we can still relax an edge, negative cycle exists
                if (dist[e.src] != INT_MAX &&                     // Source is reachable
                    dist[e.dest] > dist[e.src] + e.weight) {      // Can still relax

                    std::cout << "\n*** NEGATIVE CYCLE DETECTED ***" << std::endl;
                    std::cout << "Graph contains a negative weight cycle!" << std::endl;
                    std::cout << "Shortest path distances are not well-defined." << std::endl;
                    std::cout << "Negative cycle involves edge: (" << e.src << " -> " << e.dest
                              << ") with weight " << e.weight << std::endl;

                    // Return empty vectors to indicate negative cycle
                    return std::make_pair(std::vector<int>(), std::vector<int>());
                }
            }
        }

        return std::make_pair(dist, predecessor);                 // CHANGE: Return both
    }

    void displayShortestPathsWithPredecessors(const std::vector<int>& distances, const std::vector<int>& predecessors, int source) {
        std::cout << "Shortest paths from vertex " << source << ":" << std::endl;
        for (int i = 0; i < distances.size(); i++) {
            if (distances[i] == INT_MAX) {
                std::cout << "Vertex " << i << ": Distance = INF [p: -]" << std::endl;
            } else if (i == source) {
                std::cout << "Vertex " << i << ": Distance = " << distances[i] << " (source) [p: -]" << std::endl;
            } else {
                std::cout << "Vertex " << i << ": Distance = " << distances[i] << " [p: " << predecessors[i] << "]" << std::endl;
            }
        }
    }
    /**
     * UTILITY FUNCTION: Calculate total weight of MST
     * Used for verification and comparison purposes
     */
    int calculateTotalWeight(const vector<Edge>& mst) {
        int totalWeight = 0;                                      // Initialize sum
        for (const auto& edge : mst) {                            // For each edge in MST
            totalWeight += edge.weight;                           // Add weight to sum
        }
        return totalWeight;                                       // Return total weight
    }

    /**
     * OUTPUT FUNCTION: Print MST results for matrix implementation
     * Displays each edge and total weight
     */
    void printMSTforMatrix(const std::vector<Edge>& mst, const AdjacencyMatrix& graph) {
        std::cout << "Minimum Spanning Tree:" << std::endl;
        for (const auto& edge : mst) {                            // For each edge in MST
            std::cout << "Edge: (" << edge.src << ", " << edge.dest
                      << ") - Weight: " << edge.weight << std::endl; // Print edge details
        }
        cout << "Total weight of MST: " << calculateTotalWeight(mst) << endl; // Print total weight
    }

    /**
     * OUTPUT FUNCTION: Print MST results for list implementation
     * Displays each edge and total weight
     */
    void printMSTforList(const vector<Edge>& mst) {
        std::cout << "Minimum Spanning Tree:" << endl;
        for (const auto& edge : mst) {                            // For each edge in MST
            cout << "Edge: (" << edge.src << ", " << edge.dest
                 << ") - Weight: " << edge.weight << endl;         // Print edge details
        }
        cout << "Total weight of MST: " << calculateTotalWeight(mst) << endl; // Print total weight
    }

    /**
     * PERFORMANCE TEST FUNCTION: MST algorithms with matrix representation
     * Measures execution time and displays results
     */
    void testAdjencyMatrixMST(AdjacencyMatrix& obj, int n) {
        // PRIM'S ALGORITHM TIMING
        auto start = chrono::high_resolution_clock::now();        // Start timer
        auto mst1 = primMSTMatrix(obj);                           // Run Prim's algorithm
        auto end1 = chrono::high_resolution_clock::now();         // End timer
        cout << "Execution time of Prim's algorithm (adjacency matrix): "
             << chrono::duration_cast<chrono::microseconds>(end1 - start).count()
             << " microseconds" << endl;                          // Display timing

        // KRUSKAL'S ALGORITHM TIMING
        start = chrono::high_resolution_clock::now();             // Start timer
        auto mst2 = kruskalMSTMatrix(obj);                        // Run Kruskal's algorithm
        end1 = chrono::high_resolution_clock::now();              // End timer
        cout << "Execution time of Kruskal's algorithm (adjacency matrix): "
             << chrono::duration_cast<chrono::microseconds>(end1 - start).count()
             << " microseconds" << endl;                          // Display timing

        // DISPLAY RESULTS
        double density = 0.25;                                    // Graph density for display
        cout << "Graph size: " << n << ", Density: " << density << endl;
        cout << "Prim's algorithm (adjacency matrix):" << endl;
        printMSTforMatrix(mst1, obj);                             // Print Prim's results
        cout << "Kruskal's algorithm (adjacency matrix):" << endl;
        printMSTforMatrix(mst2, obj);                             // Print Kruskal's results
        cout << endl;
    }

    /**
     * PERFORMANCE TEST FUNCTION: SPP algorithms with matrix representation
     * Tests both Bellman-Ford and Dijkstra algorithms
     * Updated to handle negative cycle detection
     */
    void testAdjencyMatrixSPP(AdjacencyMatrix& obj, int n, int dijkstraStart) {
        // BELLMAN-FORD ALGORITHM TIMING
        auto start1 = chrono::high_resolution_clock::now();       // Start timer
        auto bellman_result = FordBellmanMatrix(obj, dijkstraStart); // Run Bellman-Ford
        auto end1 = chrono::high_resolution_clock::now();         // End timer
        cout << "Execution time of Ford-Bellman algorithm (adjacency matrix): "
             << chrono::duration_cast<chrono::microseconds>(end1 - start1).count()
             << " microseconds" << endl;                          // Display timing

        // Check if negative cycle was detected
        if (bellman_result.first.empty()) {
            cout << "Cannot compute shortest paths due to negative cycle." << endl;
            cout << "Skipping Dijkstra algorithm (not applicable with negative weights)." << endl;
            return;
        }

        // DIJKSTRA'S ALGORITHM TIMING (only if no negative cycle)
        start1 = chrono::high_resolution_clock::now();            // Start timer
        auto dijkstra_result = dijkstraMatrix(obj.getMatrix(), n, dijkstraStart); // Run Dijkstra
        end1 = chrono::high_resolution_clock::now();              // End timer
        cout << "Execution time of Dijkstra algorithm (adjacency matrix): "
             << chrono::duration_cast<chrono::microseconds>(end1 - start1).count()
             << " microseconds" << endl;                          // Display timing

        // DISPLAY BELLMAN-FORD RESULTS WITH PREDECESSORS
        cout << "Ford-Bellman algorithm (adjacency matrix):" << endl;
        displayShortestPathsWithPredecessors(bellman_result.first, bellman_result.second, dijkstraStart);
        cout << endl;

        // DISPLAY DIJKSTRA RESULTS WITH PREDECESSORS
        cout << "Dijkstra algorithm (adjacency matrix):" << endl;
        displayShortestPathsWithPredecessors(dijkstra_result.first, dijkstra_result.second, dijkstraStart);
        cout << endl;
    }
    /**
     * PERFORMANCE TEST FUNCTION: MST algorithms with adjacency list representation
     * More efficient for sparse graphs
     */
    void testAdjencyListMST(AdjacencyList& obj, int n) {
        // PRIM'S ALGORITHM TIMING
        auto start2 = chrono::high_resolution_clock::now();       // Start timer
        auto mst3 = primMSTList(obj);                             // Run Prim's algorithm
        auto end2 = chrono::high_resolution_clock::now();         // End timer
        cout << "Execution time of Prim's algorithm (adjacency list): "
             << chrono::duration_cast<chrono::microseconds>(end2 - start2).count()
             << " microseconds" << endl;                          // Display timing

        // KRUSKAL'S ALGORITHM TIMING
        start2 = chrono::high_resolution_clock::now();            // Start timer
        auto mst4 = kruskalMSTList(obj.getList(), n);             // Run Kruskal's algorithm
        end2 = chrono::high_resolution_clock::now();              // End timer
        cout << "Execution time of Kruskal's algorithm (adjacency list): "
             << chrono::duration_cast<chrono::microseconds>(end2 - start2).count()
             << " microseconds" << endl;                          // Display timing

        // DISPLAY RESULTS
        double density = 0.25;                                    // Graph density for display
        cout << "Graph size: " << n << ", Density: " << density << endl;
        cout << "Prim's algorithm (adjacency list):" << endl;
        printMSTforList(mst3);                                    // Print Prim's results
        cout << "Kruskal's algorithm (adjacency list):" << endl;
        printMSTforList(mst4);                                    // Print Kruskal's results
        cout << endl;
    }

    /**
     * PERFORMANCE TEST FUNCTION: SPP algorithms with adjacency list representation
     * More efficient for sparse graphs
     * Updated to handle negative cycle detection
     */
    void testAdjencyListSPP(AdjacencyList& obj, int n, int dijkstraStart) {
        // BELLMAN-FORD ALGORITHM TIMING
        auto start3 = chrono::high_resolution_clock::now();       // Start timer
        auto bellman_result = FordBellmanList(obj.getList(), n, dijkstraStart); // Run Bellman-Ford
        auto end3 = chrono::high_resolution_clock::now();         // End timer
        cout << "Execution time of FordBellman's algorithm (adjacency list): "
             << chrono::duration_cast<chrono::microseconds>(end3 - start3).count()
             << " microseconds" << endl;                          // Display timing

        // Check if negative cycle was detected
        if (bellman_result.first.empty()) {
            cout << "Cannot compute shortest paths due to negative cycle." << endl;
            cout << "Skipping Dijkstra algorithm (not applicable with negative weights)." << endl;
            return;
        }

        // DIJKSTRA'S ALGORITHM TIMING (only if no negative cycle)
        start3 = chrono::high_resolution_clock::now();            // Start timer
        auto dijkstra_result = dijkstraList(obj, dijkstraStart);  // Run Dijkstra
        end3 = chrono::high_resolution_clock::now();              // End timer
        cout << "Execution time of Dijkstra algorithm (adjacency list): "
             << chrono::duration_cast<chrono::microseconds>(end3 - start3).count()
             << " microseconds" << endl;                          // Display timing

        // DISPLAY BELLMAN-FORD RESULTS WITH PREDECESSORS
        cout << "FordBellman algorithm (adjacency list):" << endl;
        displayShortestPathsWithPredecessors(bellman_result.first, bellman_result.second, dijkstraStart);
        cout << endl;

        // DISPLAY DIJKSTRA RESULTS WITH PREDECESSORS
        cout << "Dijkstra algorithm (adjacency list):" << endl;
        displayShortestPathsWithPredecessors(dijkstra_result.first, dijkstra_result.second, dijkstraStart);
        cout << endl;
    }



} // namespace MST_SPP