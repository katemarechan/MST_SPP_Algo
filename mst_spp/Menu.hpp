#ifndef MENU_HPP
#define MENU_HPP

#include "MST_SPP.hpp"
#include "ReadFromFile.hpp"
#include <iostream>
#include <random>
#include <set>
#include <iomanip>
#include <limits>

using namespace std;
using namespace MST_SPP;

void mstFoo();
void sppFoo();
void dijkstraFoo();
void bellmanFordFoo();
void testAlgo();
void populateGraphs(AdjacencyMatrix& matrix, AdjacencyList& list, int vertices, int edgeCount);
void displayGraphs(const AdjacencyMatrix& matrix, const AdjacencyList& list);
void displayMSTResults(const vector<Edge>& mst_prim_matrix, const vector<Edge>& mst_prim_list,
                       const vector<Edge>& mst_kruskal_matrix, const vector<Edge>& mst_kruskal_list);
void displaySPPResults(const vector<int>& spp_dijkstra_matrix, const vector<int>& spp_dijkstra_list,
                       const vector<int>& spp_bellman_matrix, const vector<int>& spp_bellman_list,
                       int startVertex);
void displaySPPPath(const vector<int>& distances, int startVertex);

char sideMenu(string name)
{
    char choice;

    cout << "\n=== " << name << " Menu ===" << endl;
    cout << "[1] Read From File" << endl;
    cout << "[2] Generate Random Graph" << endl;
    cout << "[3] Test Algorithms (Matrix vs List)" << endl;
    cout << "[4] Display Graph Representations" << endl;
    cout << "[0] Back" << endl;
    cout << "Your Choice: ";
    cin >> choice;
    cout << endl;

    return choice;
}

void menuAlgo()
{
    int choice;

    do {
        cout << "\t=== Main Menu ===" << endl;
        cout << "[1] MST [Minimum Spanning Tree]" << endl;
        cout << "[2] SPP [Shortest Path Problem]" << endl;
        cout << "[3] Performance Comparison" << endl;
        cout << "[0] Exit" << endl;
        cout << "Please provide your choice: ";
        cin >> choice;

        switch (choice) {
            case 1:
                cout << "\n=== MST [Minimum Spanning Tree] ===" << endl;
                mstFoo();
                break;
            case 2:
                cout << "\n=== SPP [Shortest Path Problem] ===" << endl;
                sppFoo();
                break;
            case 3:
                cout << "\n=== Performance Comparison ===" << endl;
                testAlgo();
                break;
            case 0:
                cout << "Exiting program..." << endl;
                break;
            default:
                cout << "Invalid choice. Please try again." << endl;
                break;
        }
    } while (choice != 0);
}

void mstFoo() {
    char choice;
    string filename;
    AdjacencyMatrix* currentMatrix = nullptr;
    AdjacencyList* currentList = nullptr;

    do {
        choice = sideMenu("MST");

        switch (choice) {
            case '0':    //Back
                if (currentMatrix) delete currentMatrix;
                if (currentList) delete currentList;
                return;

            case '1':    //Read from file
            {
                cout << "Name of the file: ";
                cin.ignore();  // Clear input buffer
                getline(cin, filename);

                if (currentMatrix) delete currentMatrix;
                if (currentList) delete currentList;

                currentMatrix = nullptr;
                currentList = nullptr;

                bool success = FileHandler::readFromFile(filename, currentMatrix, currentList);

                if (success) {
                    cout << "\n=== Graph loaded successfully! ===" << endl;
                    cout << "Use option [4] to display graph representations" << endl;
                } else {
                    cout << "\n=== Failed to load graph! Please check the file. ===" << endl;
                    if (currentMatrix) delete currentMatrix;
                    if (currentList) delete currentList;
                    currentMatrix = nullptr;
                    currentList = nullptr;
                }
            }
                break;

            case '2':    //Generate random graph
            {
                int vertices, density_percent;
                cout << "Enter the number of vertices: ";
                cin >> vertices;
                cout << "Give the density in percentage (0-100): ";
                cin >> density_percent;

                double density = density_percent / 100.0;
                int edges = static_cast<int>(vertices * (vertices - 1) * density / 2);

                if (currentMatrix) delete currentMatrix;
                if (currentList) delete currentList;

                currentMatrix = new AdjacencyMatrix(vertices);
                currentList = new AdjacencyList(vertices);

                populateGraphs(*currentMatrix, *currentList, vertices, edges);

                cout << "\n=== Random graph generated successfully! ===" << endl;
                cout << "Graph: " << vertices << " vertices, " << edges << " edges (density: " << density * 100 << "%)" << endl;
                cout << "Use option [4] to display graph representations" << endl;
            }
                break;

            case '3':    //Test algorithms comparison
            {
                if (!currentMatrix || !currentList) {
                    cout << "Please load or generate a graph first!" << endl;
                    break;
                }

                cout << "\n=== Testing Both Representations on Current Graph ===" << endl;
                cout << "Graph has " << currentMatrix->getV() << " vertices" << endl;

                try {
                    // Test Prim's on Matrix
                    cout << "Testing Prim's algorithm on matrix..." << endl;
                    auto start = chrono::high_resolution_clock::now();
                    auto mst_prim_matrix = primMSTMatrix(*currentMatrix);
                    auto end = chrono::high_resolution_clock::now();
                    auto time_prim_matrix = chrono::duration_cast<chrono::microseconds>(end - start).count();
                    cout << "Prim's matrix completed successfully" << endl;

                    // Test Prim's on List
                    cout << "Testing Prim's algorithm on list..." << endl;
                    start = chrono::high_resolution_clock::now();
                    auto mst_prim_list = primMSTList(*currentList);
                    end = chrono::high_resolution_clock::now();
                    auto time_prim_list = chrono::duration_cast<chrono::microseconds>(end - start).count();
                    cout << "Prim's list completed successfully" << endl;

                    // Test Kruskal's on Matrix
                    cout << "Testing Kruskal's algorithm on matrix..." << endl;
                    start = chrono::high_resolution_clock::now();
                    auto mst_kruskal_matrix = kruskalMSTMatrix(*currentMatrix);
                    end = chrono::high_resolution_clock::now();
                    auto time_kruskal_matrix = chrono::duration_cast<chrono::microseconds>(end - start).count();
                    cout << "Kruskal's matrix completed successfully" << endl;

                    // Test Kruskal's on List
                    cout << "Testing Kruskal's algorithm on list..." << endl;
                    start = chrono::high_resolution_clock::now();
                    auto mst_kruskal_list = kruskalMSTList(currentList->getList(), currentMatrix->getV());
                    end = chrono::high_resolution_clock::now();
                    auto time_kruskal_list = chrono::duration_cast<chrono::microseconds>(end - start).count();
                    cout << "Kruskal's list completed successfully" << endl;

                    cout << "\n=== ALGORITHM COMPARISON RESULTS ===" << endl;
                    cout << "Prim's Algorithm:" << endl;
                    cout << "  Matrix: " << time_prim_matrix << " ?s | Total weight: " << calculateTotalWeight(mst_prim_matrix) << endl;
                    cout << "  List:   " << time_prim_list << " ?s | Total weight: " << calculateTotalWeight(mst_prim_list) << endl;

                    cout << "Kruskal's Algorithm:" << endl;
                    cout << "  Matrix: " << time_kruskal_matrix << " ?s | Total weight: " << calculateTotalWeight(mst_kruskal_matrix) << endl;
                    cout << "  List:   " << time_kruskal_list << " ?s | Total weight: " << calculateTotalWeight(mst_kruskal_list) << endl;

                    // Display detailed MST results with step-by-step edge selection
                    displayMSTResults(mst_prim_matrix, mst_prim_list, mst_kruskal_matrix, mst_kruskal_list);

                } catch (const std::exception& e) {
                    cout << "Exception caught: " << e.what() << endl;
                } catch (...) {
                    cout << "Unknown exception caught!" << endl;
                }
            }
                break;

            case '4':    //Display graph representations
            {
                if (!currentMatrix || !currentList) {
                    cout << "Please load or generate a graph first!" << endl;
                    break;
                }
                displayGraphs(*currentMatrix, *currentList);
            }
                break;

            default:
                cout << "Invalid choice. Please try again." << endl;
                break;
        }
    } while (choice != '0');

    if (currentMatrix) delete currentMatrix;
    if (currentList) delete currentList;
}

void sppFoo() {
    int choice;

    do {
        cout << "\n=== SPP [Shortest Path Problem] Algorithm Selection ===" << endl;
        cout << "[1] Dijkstra's Algorithm" << endl;
        cout << "[2] Bellman-Ford Algorithm" << endl;
        cout << "[0] Main Menu" << endl;
        cout << "Please provide your choice: ";
        cin >> choice;

        switch (choice) {
            case 1:
                cout << "\n=== Dijkstra's Algorithm ===" << endl;
                dijkstraFoo();
                break;
            case 2:
                cout << "\n=== Bellman-Ford Algorithm ===" << endl;
                bellmanFordFoo();
                break;
            case 0:
                cout << "Returning to Main Menu..." << endl;
                break;
            default:
                cout << "Invalid choice. Please try again." << endl;
                break;
        }
    } while (choice != 0);
}

void dijkstraFoo() {
    char choice;
    string filename;
    AdjacencyMatrix* currentMatrix = nullptr;
    AdjacencyList* currentList = nullptr;
    int startVertexFromFile = -1;  // Track start vertex from file

    do {
        choice = sideMenu("Dijkstra");

        switch (choice) {
            case '0':    //Back
                if (currentMatrix) delete currentMatrix;
                if (currentList) delete currentList;
                return;

            case '1':    //Read from file
            {
                cout << "Name of the file: ";
                cin.ignore();  // Clear input buffer
                getline(cin, filename);

                if (currentMatrix) delete currentMatrix;
                if (currentList) delete currentList;

                currentMatrix = nullptr;
                currentList = nullptr;
                startVertexFromFile = -1;  // Reset start vertex

                bool success = FileHandler::readFromFileForSPP(filename, currentMatrix, currentList, &startVertexFromFile);

                if (success) {
                    cout << "\n=== Graph loaded successfully! ===" << endl;
                    cout << "Use option [4] to display graph representations" << endl;
                    if (startVertexFromFile != -1) {
                        cout << "Start vertex from file: " << startVertexFromFile << endl;
                        cout << "This will be used automatically when testing algorithms." << endl;
                    } else {
                        cout << "When testing algorithms, you'll be asked for a start vertex." << endl;
                    }
                } else {
                    cout << "\n=== Failed to load graph! Please check the file. ===" << endl;
                    if (currentMatrix) delete currentMatrix;
                    if (currentList) delete currentList;
                    currentMatrix = nullptr;
                    currentList = nullptr;
                    startVertexFromFile = -1;
                }
            }
                break;

            case '2':    //Generate random graph
            {
                int vertices, density_percent;
                cout << "Enter the number of vertices: ";
                cin >> vertices;
                cout << "Give the density in percentage (0-100): ";
                cin >> density_percent;

                double density = density_percent / 100.0;
                int edges = static_cast<int>(vertices * (vertices - 1) * density / 2);

                if (currentMatrix) delete currentMatrix;
                if (currentList) delete currentList;

                currentMatrix = new AdjacencyMatrix(vertices);
                currentList = new AdjacencyList(vertices);
                startVertexFromFile = -1;  // Reset start vertex (will ask user later)

                populateGraphs(*currentMatrix, *currentList, vertices, edges);

                cout << "\n=== Random graph generated successfully! ===" << endl;
                cout << "Graph: " << vertices << " vertices, " << edges << " edges (density: " << density * 100 << "%)" << endl;
                cout << "Use option [4] to display graph representations" << endl;
                cout << "When testing algorithms, you'll be asked for a start vertex." << endl;
            }
                break;

            case '3':    //Test Dijkstra algorithms comparison
            {
                if (!currentMatrix || !currentList) {
                    cout << "Please load or generate a graph first!" << endl;
                    break;
                }

                int startVertex;
                if (startVertexFromFile != -1) {
                    // Use start vertex from file
                    startVertex = startVertexFromFile;
                    cout << "Using start vertex from file: " << startVertex << endl;
                } else {
                    // Ask user for start vertex (random graph case)
                    cout << "Enter starting vertex (0-" << currentMatrix->getV()-1 << "): ";
                    cin >> startVertex;

                    // Validate start vertex
                    if (startVertex < 0 || startVertex >= currentMatrix->getV()) {
                        cout << "Invalid start vertex! Using vertex 0 instead." << endl;
                        startVertex = 0;
                    }
                }

                cout << "\n=== Testing Dijkstra's Algorithm on Both Representations ===" << endl;
                cout << "Graph has " << currentMatrix->getV() << " vertices, starting from vertex " << startVertex << endl;

                try {
                    // Test Dijkstra on Matrix
                    cout << "Testing Dijkstra's algorithm on matrix..." << endl;
                    auto start = chrono::high_resolution_clock::now();
                    auto dijkstra_matrix_result = dijkstraMatrix(currentMatrix->getMatrix(), currentMatrix->getV(), startVertex);
                    auto end = chrono::high_resolution_clock::now();
                    auto time_dijkstra_matrix = chrono::duration_cast<chrono::microseconds>(end - start).count();
                    cout << "Dijkstra's matrix completed successfully" << endl;

                    // Test Dijkstra on List
                    cout << "Testing Dijkstra's algorithm on list..." << endl;
                    start = chrono::high_resolution_clock::now();
                    auto dijkstra_list_result = dijkstraList(*currentList, startVertex);
                    end = chrono::high_resolution_clock::now();
                    auto time_dijkstra_list = chrono::duration_cast<chrono::microseconds>(end - start).count();
                    cout << "Dijkstra's list completed successfully" << endl;

                    // Extract distances from the pairs
                    vector<int> spp_dijkstra_matrix = dijkstra_matrix_result.first;
                    vector<int> spp_dijkstra_list = dijkstra_list_result.first;

                    cout << "\n=== DIJKSTRA'S ALGORITHM COMPARISON RESULTS ===" << endl;
                    cout << "Matrix: " << time_dijkstra_matrix << " ?s" << endl;
                    cout << "List:   " << time_dijkstra_list << " ?s" << endl;

                    // Display Dijkstra results with predecessors
                    cout << "\n--- Dijkstra's Algorithm (Matrix) ---" << endl;
                    displayShortestPathsWithPredecessors(dijkstra_matrix_result.first, dijkstra_matrix_result.second, startVertex);

                    cout << "\n--- Dijkstra's Algorithm (List) ---" << endl;
                    displayShortestPathsWithPredecessors(dijkstra_list_result.first, dijkstra_list_result.second, startVertex);

                    // Compare results
                    bool dijkstra_same = true;
                    if (spp_dijkstra_matrix.size() != spp_dijkstra_list.size()) {
                        dijkstra_same = false;
                    } else {
                        for (int i = 0; i < spp_dijkstra_matrix.size(); i++) {
                            if (spp_dijkstra_matrix[i] != spp_dijkstra_list[i]) {
                                dijkstra_same = false;
                                break;
                            }
                        }
                    }

                    cout << "\n=== ALGORITHM VERIFICATION ===" << endl;
                    cout << "Dijkstra Matrix vs List: " << (dijkstra_same ? "? Results match" : "? Results differ") << endl;

                } catch (const std::exception& e) {
                    cout << "Exception caught: " << e.what() << endl;
                } catch (...) {
                    cout << "Unknown exception caught!" << endl;
                }
            }
                break;

            case '4':    //Display graph representations
            {
                if (!currentMatrix || !currentList) {
                    cout << "Please load or generate a graph first!" << endl;
                    break;
                }
                displayGraphs(*currentMatrix, *currentList);
                if (startVertexFromFile != -1) {
                    cout << "\nStart vertex from file: " << startVertexFromFile << endl;
                }
            }
                break;

            default:
                cout << "Invalid choice. Please try again." << endl;
                break;
        }
    } while (choice != '0');

    if (currentMatrix) delete currentMatrix;
    if (currentList) delete currentList;
}

void bellmanFordFoo() {
    char choice;
    string filename;
    AdjacencyMatrix* currentMatrix = nullptr;
    AdjacencyList* currentList = nullptr;
    int startVertexFromFile = -1;  // Track start vertex from file

    do {
        choice = sideMenu("Bellman-Ford");

        switch (choice) {
            case '0':    //Back
                if (currentMatrix) delete currentMatrix;
                if (currentList) delete currentList;
                return;

            case '1':    //Read from file
            {
                cout << "Name of the file: ";
                cin.ignore();  // Clear input buffer
                getline(cin, filename);

                if (currentMatrix) delete currentMatrix;
                if (currentList) delete currentList;

                currentMatrix = nullptr;
                currentList = nullptr;
                startVertexFromFile = -1;  // Reset start vertex

                bool success = FileHandler::readFromFileForSPP(filename, currentMatrix, currentList, &startVertexFromFile);

                if (success) {
                    cout << "\n=== Graph loaded successfully! ===" << endl;
                    cout << "Use option [4] to display graph representations" << endl;
                    if (startVertexFromFile != -1) {
                        cout << "Start vertex from file: " << startVertexFromFile << endl;
                        cout << "This will be used automatically when testing algorithms." << endl;
                    } else {
                        cout << "When testing algorithms, you'll be asked for a start vertex." << endl;
                    }
                } else {
                    cout << "\n=== Failed to load graph! Please check the file. ===" << endl;
                    if (currentMatrix) delete currentMatrix;
                    if (currentList) delete currentList;
                    currentMatrix = nullptr;
                    currentList = nullptr;
                    startVertexFromFile = -1;
                }
            }
                break;

            case '2':    //Generate random graph
            {
                int vertices, density_percent;
                cout << "Enter the number of vertices: ";
                cin >> vertices;
                cout << "Give the density in percentage (0-100): ";
                cin >> density_percent;

                double density = density_percent / 100.0;
                int edges = static_cast<int>(vertices * (vertices - 1) * density / 2);

                if (currentMatrix) delete currentMatrix;
                if (currentList) delete currentList;

                currentMatrix = new AdjacencyMatrix(vertices);
                currentList = new AdjacencyList(vertices);
                startVertexFromFile = -1;  // Reset start vertex (will ask user later)

                populateGraphs(*currentMatrix, *currentList, vertices, edges);

                cout << "\n=== Random graph generated successfully! ===" << endl;
                cout << "Graph: " << vertices << " vertices, " << edges << " edges (density: " << density * 100 << "%)" << endl;
                cout << "Use option [4] to display graph representations" << endl;
                cout << "When testing algorithms, you'll be asked for a start vertex." << endl;
            }
                break;

            case '3':    //Test Bellman-Ford algorithms comparison
            {
                if (!currentMatrix || !currentList) {
                    cout << "Please load or generate a graph first!" << endl;
                    break;
                }

                int startVertex;
                if (startVertexFromFile != -1) {
                    // Use start vertex from file
                    startVertex = startVertexFromFile;
                    cout << "Using start vertex from file: " << startVertex << endl;
                } else {
                    // Ask user for start vertex (random graph case)
                    cout << "Enter starting vertex (0-" << currentMatrix->getV()-1 << "): ";
                    cin >> startVertex;

                    // Validate start vertex
                    if (startVertex < 0 || startVertex >= currentMatrix->getV()) {
                        cout << "Invalid start vertex! Using vertex 0 instead." << endl;
                        startVertex = 0;
                    }
                }

                cout << "\n=== Testing Bellman-Ford Algorithm on Both Representations ===" << endl;
                cout << "Graph has " << currentMatrix->getV() << " vertices, starting from vertex " << startVertex << endl;

                try {
                    // Test Bellman-Ford on Matrix
                    cout << "Testing Bellman-Ford algorithm on matrix..." << endl;
                    auto start = chrono::high_resolution_clock::now();
                    auto bellman_matrix_result = FordBellmanMatrix(*currentMatrix, startVertex);
                    auto end = chrono::high_resolution_clock::now();
                    auto time_bellman_matrix = chrono::duration_cast<chrono::microseconds>(end - start).count();

                    // Check if negative cycle detected in matrix implementation
                    if (bellman_matrix_result.first.empty()) {
                        cout << "Bellman-Ford matrix: NEGATIVE CYCLE DETECTED - algorithm terminated" << endl;
                        cout << "Skipping list implementation due to negative cycle." << endl;
                        cout << "Cannot compute shortest paths when negative cycles exist." << endl;
                        break;  // Exit this test case
                    } else {
                        cout << "Bellman-Ford matrix completed successfully" << endl;
                    }

                    // Test Bellman-Ford on List (only if matrix didn't find negative cycle)
                    cout << "Testing Bellman-Ford algorithm on list..." << endl;
                    start = chrono::high_resolution_clock::now();
                    auto bellman_list_result = FordBellmanList(currentList->getList(), currentList->getV(), startVertex);
                    end = chrono::high_resolution_clock::now();
                    auto time_bellman_list = chrono::duration_cast<chrono::microseconds>(end - start).count();

                    // Check if negative cycle detected in list implementation
                    if (bellman_list_result.first.empty()) {
                        cout << "Bellman-Ford list: NEGATIVE CYCLE DETECTED - algorithm terminated" << endl;
                        cout << "Cannot compute shortest paths when negative cycles exist." << endl;
                        break;  // Exit this test case
                    } else {
                        cout << "Bellman-Ford list completed successfully" << endl;
                    }

                    cout << "\n=== BELLMAN-FORD ALGORITHM COMPARISON RESULTS ===" << endl;
                    cout << "Matrix: " << time_bellman_matrix << " ?s" << endl;
                    cout << "List:   " << time_bellman_list << " ?s" << endl;

                    // Display Bellman-Ford results with predecessors (only if both succeeded)
                    cout << "\n--- Bellman-Ford Algorithm (Matrix) ---" << endl;
                    displayShortestPathsWithPredecessors(bellman_matrix_result.first, bellman_matrix_result.second, startVertex);

                    cout << "\n--- Bellman-Ford Algorithm (List) ---" << endl;
                    displayShortestPathsWithPredecessors(bellman_list_result.first, bellman_list_result.second, startVertex);

                    // Compare results (only if both succeeded)
                    bool bellman_same = true;
                    if (bellman_matrix_result.first.size() != bellman_list_result.first.size()) {
                        bellman_same = false;
                    } else {
                        for (int i = 0; i < bellman_matrix_result.first.size(); i++) {
                            if (bellman_matrix_result.first[i] != bellman_list_result.first[i]) {
                                bellman_same = false;
                                break;
                            }
                        }
                    }

                    cout << "\n=== ALGORITHM VERIFICATION ===" << endl;
                    cout << "Bellman-Ford Matrix vs List: " << (bellman_same ? "? Results match" : "? Results differ") << endl;

                } catch (const std::exception& e) {
                    cout << "Exception caught: " << e.what() << endl;
                } catch (...) {
                    cout << "Unknown exception caught!" << endl;
                }
            }
                break;

            case '4':    //Display graph representations
            {
                if (!currentMatrix || !currentList) {
                    cout << "Please load or generate a graph first!" << endl;
                    break;
                }
                displayGraphs(*currentMatrix, *currentList);
                if (startVertexFromFile != -1) {
                    cout << "\nStart vertex from file: " << startVertexFromFile << endl;
                }
            }
                break;

            default:
                cout << "Invalid choice. Please try again." << endl;
                break;
        }
    } while (choice != '0');

    if (currentMatrix) delete currentMatrix;
    if (currentList) delete currentList;
}

void testAlgo() {
    cout << "\n=== Complete Performance Analysis ===" << endl;
    cout << "Testing all algorithms on the same graphs with different representations\n" << endl;

    vector<int> sizes = {50, 100, 200};
    vector<double> densities = {0.25, 0.5, 0.99};

    cout << "Graph Size\tDensity\tAlgorithm\t\tMatrix (?s)\tList (?s)\tRatio (M/L)" << endl;
    cout << "????????????????????????????????????????????????????????????????????????" << endl;

    for (int size : sizes) {
        for (double density : densities) {
            int edges = static_cast<int>(size * (size - 1) * density / 2);

            AdjacencyMatrix mat(size);
            AdjacencyList list(size);
            populateGraphs(mat, list, size, edges);

            // MST ALGORITHMS
            // Test Prim's Algorithm
            auto start = chrono::high_resolution_clock::now();
            auto prim_matrix_result = primMSTMatrix(mat);
            auto end = chrono::high_resolution_clock::now();
            auto time_prim_matrix = chrono::duration_cast<chrono::microseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            auto prim_list_result = primMSTList(list);
            end = chrono::high_resolution_clock::now();
            auto time_prim_list = chrono::duration_cast<chrono::microseconds>(end - start).count();

            // Test Kruskal's Algorithm
            start = chrono::high_resolution_clock::now();
            auto kruskal_matrix_result = kruskalMSTMatrix(mat);
            end = chrono::high_resolution_clock::now();
            auto time_kruskal_matrix = chrono::duration_cast<chrono::microseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            auto kruskal_list_result = kruskalMSTList(list.getList(), size);
            end = chrono::high_resolution_clock::now();
            auto time_kruskal_list = chrono::duration_cast<chrono::microseconds>(end - start).count();

            // SPP ALGORITHMS
            int startVertex = 0;

            // Test Dijkstra's Algorithm
            start = chrono::high_resolution_clock::now();
            auto dijkstra_matrix_pair = dijkstraMatrix(mat.getMatrix(), size, startVertex);
            end = chrono::high_resolution_clock::now();
            auto time_dijkstra_matrix = chrono::duration_cast<chrono::microseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            auto dijkstra_list_pair = dijkstraList(list, startVertex);
            end = chrono::high_resolution_clock::now();
            auto time_dijkstra_list = chrono::duration_cast<chrono::microseconds>(end - start).count();

            // Extract distances from Dijkstra pairs
            vector<int> dijkstra_matrix_result = dijkstra_matrix_pair.first;
            vector<int> dijkstra_list_result = dijkstra_list_pair.first;

            // Test Bellman-Ford Algorithm with negative cycle detection
            start = chrono::high_resolution_clock::now();
            auto bellman_matrix_pair = FordBellmanMatrix(mat, startVertex);
            end = chrono::high_resolution_clock::now();
            auto time_bellman_matrix = chrono::duration_cast<chrono::microseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            auto bellman_list_pair = FordBellmanList(list.getList(), size, startVertex);
            end = chrono::high_resolution_clock::now();
            auto time_bellman_list = chrono::duration_cast<chrono::microseconds>(end - start).count();

            // Check for negative cycles in Bellman-Ford results
            bool bellman_matrix_has_neg_cycle = bellman_matrix_pair.first.empty();
            bool bellman_list_has_neg_cycle = bellman_list_pair.first.empty();
            bool any_negative_cycle = bellman_matrix_has_neg_cycle || bellman_list_has_neg_cycle;

            // Extract distances from Bellman-Ford pairs (if no negative cycle)
            vector<int> bellman_matrix_result, bellman_list_result;
            if (!bellman_matrix_has_neg_cycle) {
                bellman_matrix_result = bellman_matrix_pair.first;
            }
            if (!bellman_list_has_neg_cycle) {
                bellman_list_result = bellman_list_pair.first;
            }

            // Calculate performance ratios
            double ratio_prim = (time_prim_list > 0) ? (double)time_prim_matrix / time_prim_list : 0;
            double ratio_kruskal = (time_kruskal_list > 0) ? (double)time_kruskal_matrix / time_kruskal_list : 0;
            double ratio_dijkstra = (time_dijkstra_list > 0) ? (double)time_dijkstra_matrix / time_dijkstra_list : 0;

            // OUTPUT RESULTS
            cout << size << "\t\t" << density << "\tPrim's MST\t\t" << time_prim_matrix << "\t\t" << time_prim_list << "\t\t" << fixed << setprecision(2) << ratio_prim << endl;
            cout << size << "\t\t" << density << "\tKruskal's MST\t\t" << time_kruskal_matrix << "\t\t" << time_kruskal_list << "\t\t" << ratio_kruskal << endl;
            cout << size << "\t\t" << density << "\tDijkstra's SPP\t\t" << time_dijkstra_matrix << "\t\t" << time_dijkstra_list << "\t\t" << ratio_dijkstra << endl;

            // Handle Bellman-Ford results (check for negative cycles)
            if (any_negative_cycle) {
                cout << size << "\t\t" << density << "\tBellman-Ford SPP\t";

                if (bellman_matrix_has_neg_cycle) {
                    cout << "NEG_CYCLE\t\t";
                } else {
                    cout << time_bellman_matrix << "\t\t";
                }

                if (bellman_list_has_neg_cycle) {
                    cout << "NEG_CYCLE\t\t";
                } else {
                    cout << time_bellman_list << "\t\t";
                }

                cout << "N/A" << endl;

                // Display negative cycle warning for this graph configuration
                cout << "\t\t\t*** NEGATIVE CYCLE DETECTED in " << size << "V/" << density << "D graph ***" << endl;
            } else {
                // Both Bellman-Ford implementations succeeded
                double ratio_bellman = (time_bellman_list > 0) ? (double)time_bellman_matrix / time_bellman_list : 0;
                cout << size << "\t\t" << density << "\tBellman-Ford SPP\t" << time_bellman_matrix << "\t\t" << time_bellman_list << "\t\t" << ratio_bellman << endl;
            }

            cout << "????????????????????????????????????????????????????????????????????????" << endl;

            // ALGORITHM VERIFICATION for this graph size/density
            cout << "\t\t\t=== Verification for " << size << "V/" << density << "D ===" << endl;

            // Verify MST results
            bool prim_same = (calculateTotalWeight(prim_matrix_result) == calculateTotalWeight(prim_list_result));
            bool kruskal_same = (calculateTotalWeight(kruskal_matrix_result) == calculateTotalWeight(kruskal_list_result));

            cout << "\t\t\tPrim's Matrix vs List: " << (prim_same ? "? Weights match" : "? Weights differ") << endl;
            cout << "\t\t\tKruskal's Matrix vs List: " << (kruskal_same ? "? Weights match" : "? Weights differ") << endl;

            // Verify Dijkstra results
            bool dijkstra_same = true;
            if (dijkstra_matrix_result.size() != dijkstra_list_result.size()) {
                dijkstra_same = false;
            } else {
                for (int i = 0; i < dijkstra_matrix_result.size(); i++) {
                    if (dijkstra_matrix_result[i] != dijkstra_list_result[i]) {
                        dijkstra_same = false;
                        break;
                    }
                }
            }
            cout << "\t\t\tDijkstra Matrix vs List: " << (dijkstra_same ? "? Results match" : "? Results differ") << endl;

            // Verify Bellman-Ford results (only if both succeeded)
            if (!any_negative_cycle) {
                bool bellman_same = true;
                if (bellman_matrix_result.size() != bellman_list_result.size()) {
                    bellman_same = false;
                } else {
                    for (int i = 0; i < bellman_matrix_result.size(); i++) {
                        if (bellman_matrix_result[i] != bellman_list_result[i]) {
                            bellman_same = false;
                            break;
                        }
                    }
                }
                cout << "\t\t\tBellman-Ford Matrix vs List: " << (bellman_same ? "? Results match" : "? Results differ") << endl;
            } else {
                cout << "\t\t\tBellman-Ford Matrix vs List: Cannot compare (negative cycle)" << endl;
            }

            cout << "????????????????????????????????????????????????????????????????????????" << endl;
        }

        cout << "\n"; // Extra spacing between different graph sizes
    }
}

// Display functions for required output formats
void displayGraphs(const AdjacencyMatrix& adjMatrix, const AdjacencyList& adjList) {
    cout << "\n=== GRAPH REPRESENTATIONS ===" << endl;

    cout << "\n--- Adjacency Matrix ---" << endl;
    cout << "   ";
    for (int i = 0; i < adjMatrix.getV(); i++) {
        cout << setw(4) << i;
    }
    cout << endl;

    for (int i = 0; i < adjMatrix.getV(); i++) {
        cout << setw(2) << i << " ";
        for (int j = 0; j < adjMatrix.getV(); j++) {
            cout << setw(4) << adjMatrix.getMatrix()[i][j];
        }
        cout << endl;
    }

    cout << "\n--- Adjacency List ---" << endl;
    for (int i = 0; i < adjList.getV(); i++) {
        cout << "Vertex " << i << ": ";
        for (const auto& edge : adjList.getList()[i]) {
            cout << "(" << edge.dest << "," << edge.weight << ") ";
        }
        cout << endl;
    }
}

// Display function to show step-by-step MST construction
void displayMSTResults(const vector<Edge>& mst_prim_matrix, const vector<Edge>& mst_prim_list,
                       const vector<Edge>& mst_kruskal_matrix, const vector<Edge>& mst_kruskal_list) {
    cout << "\n=== MST RESULTS DISPLAY ===" << endl;

    cout << "\n--- Prim's Algorithm (Matrix) ---" << endl;
    cout << "Edges chosen step by step:" << endl;
    for (int i = 0; i < mst_prim_matrix.size(); i++) {
        const auto& edge = mst_prim_matrix[i];
        cout << "Step " << (i + 1) << ": Edge (" << edge.src << ", " << edge.dest << ") - Weight: " << edge.weight << endl;
    }
    cout << "Total Weight: " << calculateTotalWeight(mst_prim_matrix) << endl;

    cout << "\n--- Prim's Algorithm (List) ---" << endl;
    cout << "Edges chosen step by step:" << endl;
    for (int i = 0; i < mst_prim_list.size(); i++) {
        const auto& edge = mst_prim_list[i];
        cout << "Step " << (i + 1) << ": Edge (" << edge.src << ", " << edge.dest << ") - Weight: " << edge.weight << endl;
    }
    cout << "Total Weight: " << calculateTotalWeight(mst_prim_list) << endl;

    cout << "\n--- Kruskal's Algorithm (Matrix) ---" << endl;
    cout << "Edges chosen step by step:" << endl;
    for (int i = 0; i < mst_kruskal_matrix.size(); i++) {
        const auto& edge = mst_kruskal_matrix[i];
        cout << "Step " << (i + 1) << ": Edge (" << edge.src << ", " << edge.dest << ") - Weight: " << edge.weight << endl;
    }
    cout << "Total Weight: " << calculateTotalWeight(mst_kruskal_matrix) << endl;

    cout << "\n--- Kruskal's Algorithm (List) ---" << endl;
    cout << "Edges chosen step by step:" << endl;
    for (int i = 0; i < mst_kruskal_list.size(); i++) {
        const auto& edge = mst_kruskal_list[i];
        cout << "Step " << (i + 1) << ": Edge (" << edge.src << ", " << edge.dest << ") - Weight: " << edge.weight << endl;
    }
    cout << "Total Weight: " << calculateTotalWeight(mst_kruskal_list) << endl;

    // Compare Prim's implementations
    bool prim_same = true;
    if (mst_prim_matrix.size() != mst_prim_list.size()) {
        prim_same = false;
    } else {
        for (int i = 0; i < mst_prim_matrix.size(); i++) {
            if (mst_prim_matrix[i].src != mst_prim_list[i].src ||
                mst_prim_matrix[i].dest != mst_prim_list[i].dest ||
                mst_prim_matrix[i].weight != mst_prim_list[i].weight) {
                prim_same = false;
                break;
            }
        }
    }

    // Compare Kruskal's implementations
    bool kruskal_same = true;
    if (mst_kruskal_matrix.size() != mst_kruskal_list.size()) {
        kruskal_same = false;
    } else {
        for (int i = 0; i < mst_kruskal_matrix.size(); i++) {
            if (mst_kruskal_matrix[i].src != mst_kruskal_list[i].src ||
                mst_kruskal_matrix[i].dest != mst_kruskal_list[i].dest ||
                mst_kruskal_matrix[i].weight != mst_kruskal_list[i].weight) {
                kruskal_same = false;
                break;
            }
        }
    }

    cout << "\n=== ALGORITHM VERIFICATION ===" << endl;
    cout << "Prim's Matrix vs List: " << (prim_same ? "? Results match" : "? Results differ") << endl;
    cout << "Kruskal's Matrix vs List: " << (kruskal_same ? "? Results match" : "? Results differ") << endl;
}

void displaySPPResults(const vector<int>& spp_dijkstra_matrix, const vector<int>& spp_dijkstra_list,
                       const vector<int>& spp_bellman_matrix, const vector<int>& spp_bellman_list,
                       int startVertex) {
    cout << "\n=== SPP RESULTS DISPLAY ===" << endl;

    cout << "\n--- Dijkstra's Algorithm (Matrix) ---" << endl;
    cout << "Shortest paths from vertex " << startVertex << ":" << endl;
    for (int i = 0; i < spp_dijkstra_matrix.size(); i++) {
        if (i == startVertex) {
            cout << "Vertex " << i << ": Distance = 0 (source)" << endl;
        } else if (spp_dijkstra_matrix[i] == INT_MAX) {
            cout << "Vertex " << i << ": Distance = ? (unreachable)" << endl;
        } else {
            cout << "Vertex " << i << ": Distance = " << spp_dijkstra_matrix[i] << endl;
        }
    }

    cout << "\n--- Dijkstra's Algorithm (List) ---" << endl;
    cout << "Shortest paths from vertex " << startVertex << ":" << endl;
    for (int i = 0; i < spp_dijkstra_list.size(); i++) {
        if (i == startVertex) {
            cout << "Vertex " << i << ": Distance = 0 (source)" << endl;
        } else if (spp_dijkstra_list[i] == INT_MAX) {
            cout << "Vertex " << i << ": Distance = ? (unreachable)" << endl;
        } else {
            cout << "Vertex " << i << ": Distance = " << spp_dijkstra_list[i] << endl;
        }
    }

    cout << "\n--- Bellman-Ford Algorithm (Matrix) ---" << endl;
    cout << "Shortest paths from vertex " << startVertex << ":" << endl;
    for (int i = 0; i < spp_bellman_matrix.size(); i++) {
        if (i == startVertex) {
            cout << "Vertex " << i << ": Distance = 0 (source)" << endl;
        } else if (spp_bellman_matrix[i] == INT_MAX) {
            cout << "Vertex " << i << ": Distance = ? (unreachable)" << endl;
        } else {
            cout << "Vertex " << i << ": Distance = " << spp_bellman_matrix[i] << endl;
        }
    }

    cout << "\n--- Bellman-Ford Algorithm (List) ---" << endl;
    cout << "Shortest paths from vertex " << startVertex << ":" << endl;
    for (int i = 0; i < spp_bellman_list.size(); i++) {
        if (i == startVertex) {
            cout << "Vertex " << i << ": Distance = 0 (source)" << endl;
        } else if (spp_bellman_list[i] == INT_MAX) {
            cout << "Vertex " << i << ": Distance = ? (unreachable)" << endl;
        } else {
            cout << "Vertex " << i << ": Distance = " << spp_bellman_list[i] << endl;
        }
    }

    // Comparison section (like your MST comparison)
    cout << "\n=== ALGORITHM VERIFICATION ===" << endl;

    // Compare Dijkstra implementations
    bool dijkstra_same = true;
    if (spp_dijkstra_matrix.size() != spp_dijkstra_list.size()) {
        dijkstra_same = false;
    } else {
        for (int i = 0; i < spp_dijkstra_matrix.size(); i++) {
            if (spp_dijkstra_matrix[i] != spp_dijkstra_list[i]) {
                dijkstra_same = false;
                break;
            }
        }
    }

    // Compare Bellman-Ford implementations
    bool bellman_same = true;
    if (spp_bellman_matrix.size() != spp_bellman_list.size()) {
        bellman_same = false;
    } else {
        for (int i = 0; i < spp_bellman_matrix.size(); i++) {
            if (spp_bellman_matrix[i] != spp_bellman_list[i]) {
                bellman_same = false;
                break;
            }
        }
    }

    cout << "Dijkstra Matrix vs List: " << (dijkstra_same ? "? Results match" : "? Results differ") << endl;
    cout << "Bellman-Ford Matrix vs List: " << (bellman_same ? "? Results match" : "? Results differ") << endl;
}

void displaySPPPath(const vector<int>& distances, int startVertex) {
    cout << "\n--- Path Summary ---" << endl;
    cout << "Paths from vertex " << startVertex << " in sequence:" << endl;
    for (int i = 0; i < distances.size(); i++) {
        if (i != startVertex) {
            cout << startVertex << " ? " << i << ": ";
            if (distances[i] == INT_MAX) {
                cout << "No path" << endl;
            } else {
                cout << "Distance " << distances[i] << endl;
            }
        }
    }
}

// Helper function to populate both representations with the same graph
void populateGraphs(AdjacencyMatrix& matrix, AdjacencyList& list, int vertices, int edgeCount) {
    random_device rd;
    mt19937 gen(rd());

    // Generate the same edges for both representations
    vector<tuple<int, int, int>> edges;
    set<pair<int, int>> edgeSet;

    // Generate unique edges
    while (edges.size() < edgeCount) {
        int src = gen() % vertices;
        int dest = gen() % vertices;

        if (src != dest) {
            pair<int, int> edge = {min(src, dest), max(src, dest)};
            if (edgeSet.find(edge) == edgeSet.end()) {
                int weight = gen() % 100 + 1;
                edges.push_back({src, dest, weight});
                edgeSet.insert(edge);
            }
        }
    }

    // Add the same edges to both representations
    for (const auto& edge : edges) {
        int src, dest, weight;
        tie(src, dest, weight) = edge;
        matrix.addEdge(src, dest, weight);
        list.addEdge(src, dest, weight);
    }
}

#endif // MENU_HPP