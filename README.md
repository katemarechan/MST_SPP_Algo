# MST_SPP_Algo
MST, SPP algorithms

This project contains the next algorithms and functionalities:

Reads the graphs, directed and not from the files, generates a graphs with the parameters that the user wants (number or E, V and density)
Performs tests for the next algorithms: Prim, Kruskal, Dijkstra, Bellman-Ford

The specifics for each algorithm:
Prim: works on min-heap with the priority queue
Kruskal: Has Union by Rank implemented
Dijkstra: also works on a min-heap for better perfomance
Bellman-Ford: Has a negative cycle detection, but doesnt stop as it detects the first cycle itself, the algorithm works V times and detects on the Vth iteration.
If the cycle is detected on the matrix -- it will automatically not work for the list.

All algorithms work on lists and matrixes.
For MST algorithms it shows the end results and the weights of the edges, while SPP shows the predecessors next to each edge
Theres an option called 'tests' where you can test the time it takes for the algorithms to work, for densities 25%, 50% and 99%

The program also shows how the matrix and the list looks for each graph created if asked
