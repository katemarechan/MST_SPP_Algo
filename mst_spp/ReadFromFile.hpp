#ifndef READFROMFILE_HPP
#define READFROMFILE_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "MST_SPP.hpp"

namespace MST_SPP {

    class FileHandler {
    public:

        // Main method that detects file format and creates appropriate graph
        static bool readFromFile(const std::string &filename, AdjacencyMatrix*& matrix, AdjacencyList*& list, bool forSPP = false, int* startVertex = nullptr) {

            std::vector<std::string> pathsToTry = {
                    filename,
                    filename + ".txt",
                    "./" + filename,
                    "./" + filename + ".txt"
            };

            std::ifstream infile;
            std::string actualPath;

            for (const auto& path : pathsToTry) {
                infile.open(path);
                if (infile.is_open()) {
                    actualPath = path;
                    std::cout << "Found file: " << actualPath << std::endl;
                    break;
                }
                infile.clear();
            }

            if (!infile.is_open()) {
                std::cerr << "Error: Could not find file. Tried:" << std::endl;
                for (const auto& path : pathsToTry) {
                    std::cerr << "  - " << path << std::endl;
                }
                return false;
            }

            // Try to read the first line - could be 2 or 3 numbers
            std::string firstLine;
            std::getline(infile, firstLine);
            std::istringstream iss(firstLine);

            int edges, vertices, startVertexFromFile = -1;
            iss >> edges >> vertices;

            // Check if there's a third number (start vertex for SPP)
            if (iss >> startVertexFromFile) {
                std::cout << "Detected SPP format: " << edges << " edges, " << vertices << " vertices, start vertex: " << startVertexFromFile << std::endl;
                if (startVertex) *startVertex = startVertexFromFile;  // Return start vertex if pointer provided
            } else {
                std::cout << "Detected MST format: " << edges << " edges, " << vertices << " vertices" << std::endl;
                if (startVertex) *startVertex = -1;  // No start vertex for MST
            }

            // Validate input
            if (vertices <= 0 || edges < 0) {
                std::cerr << "Error: Invalid graph parameters. Vertices: " << vertices << ", Edges: " << edges << std::endl;
                infile.close();
                return false;
            }

            std::cout << "Creating graph with " << vertices << " vertices..." << std::endl;

            // Create new objects with the correct size
            matrix = new AdjacencyMatrix(vertices);
            list = new AdjacencyList(vertices);

            int src, dest, weight;
            int edgeCount = 0;

            while (infile >> src >> dest >> weight) {
                // Validate vertex indices
                if (src < 0 || src >= vertices || dest < 0 || dest >= vertices) {
                    std::cerr << "Warning: Invalid edge (" << src << "," << dest << ") - vertex indices out of range [0," << vertices-1 << "] - skipping" << std::endl;
                    continue;
                }

                if (forSPP || startVertexFromFile != -1) {
                    // For shortest path algorithms, use directed edges
                    matrix->getMatrix()[src][dest] = weight;
                    list->getList()[src].push_back({src, dest, weight});
                } else {
                    // For MST algorithms, use undirected edges
                    matrix->addEdge(src, dest, weight);
                    list->addEdge(src, dest, weight);
                }

                edgeCount++;
            }

            infile.close();

            if (startVertexFromFile != -1) {
                std::cout << "Start vertex for shortest path algorithms: " << startVertexFromFile << std::endl;
            }

            std::cout << "Successfully loaded graph with " << vertices << " vertices and "
                      << edgeCount << " edges." << std::endl;
            return true;
        }

        // Convenience methods for specific algorithm types
        static bool readFromFileForMST(const std::string &filename, AdjacencyMatrix*& matrix, AdjacencyList*& list) {
            return readFromFile(filename, matrix, list, false, nullptr);
        }

        static bool readFromFileForSPP(const std::string &filename, AdjacencyMatrix*& matrix, AdjacencyList*& list, int* startVertex = nullptr) {
            return readFromFile(filename, matrix, list, true, startVertex);
        }

        // Alternative method for reading undirected graphs (MST)
        static bool readFromFileUndirected(const std::string &filename, AdjacencyMatrix*& matrix, AdjacencyList*& list) {

            std::vector<std::string> pathsToTry = {
                    filename,
                    filename + ".txt",
                    "./" + filename,
                    "./" + filename + ".txt"
            };

            std::ifstream infile;
            std::string actualPath;

            for (const auto& path : pathsToTry) {
                infile.open(path);
                if (infile.is_open()) {
                    actualPath = path;
                    std::cout << "Found file: " << actualPath << std::endl;
                    break;
                }
                infile.clear();
            }

            if (!infile.is_open()) {
                std::cerr << "Error: Could not find file. Tried:" << std::endl;
                for (const auto& path : pathsToTry) {
                    std::cerr << "  - " << path << std::endl;
                }
                return false;
            }

            int edges, vertices;
            if (!(infile >> edges >> vertices)) {
                std::cerr << "Error: Could not read graph parameters from file" << std::endl;
                infile.close();
                return false;
            }

            // Validate input
            if (vertices <= 0 || edges < 0) {
                std::cerr << "Error: Invalid graph parameters. Vertices: " << vertices << ", Edges: " << edges << std::endl;
                infile.close();
                return false;
            }

            std::cout << "Creating graph with " << vertices << " vertices..." << std::endl;

            // Create new objects with the correct size
            matrix = new AdjacencyMatrix(vertices);
            list = new AdjacencyList(vertices);

            int src, dest, weight;
            int edgeCount = 0;

            while (infile >> src >> dest >> weight) {
                // Validate vertex indices
                if (src < 0 || src >= vertices || dest < 0 || dest >= vertices) {
                    std::cerr << "Warning: Invalid edge (" << src << "," << dest << ") - vertex indices out of range [0," << vertices-1 << "] - skipping" << std::endl;
                    continue;
                }

                // For MST, make edges undirected
                matrix->addEdge(src, dest, weight);
                list->addEdge(src, dest, weight);
                edgeCount++;
            }

            infile.close();

            std::cout << "Successfully loaded graph with " << vertices << " vertices and "
                      << edgeCount << " edges." << std::endl;
            return true;
        }
    };

} // namespace MST_SPP

#endif // READFROMFILE_HPP