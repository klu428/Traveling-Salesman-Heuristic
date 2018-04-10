# Traveling Salesman Heuristic

This program is a heuristic for solving the Traveling Salesman Problem. For explanation please read about the problem on Wikipedia: https://en.wikipedia.org/wiki/Traveling_salesman_problem.

## Christofides Algorithm

For data sets of all sizes we first use the Christofides algorithm.
This is a heuristic algorithm that first finds a minimum spanning tree of a given graph.
Then, it looks just at the set of vertices with odd degrees. Using those vertices, it finds a
minimum-weight “perfect matching”, which is the set of edges without common vertices.
It then combines the edges from the minimum spanning tree and the edges of the
perfect matching to create a connected multigraph. In this multigraph, every vertex has
an even degree. Then, it finds an Euler circuit in that multigraph. Finally, it uses that
circuit to find a Hamiltonian circuit by skipping any repeated vertices.

## 2-opt Algorithm

For smaller data sets (under 1000), we use the 2-Opt algorithm to improve our solution.
This is an optimization algorithm that is applied to an already generated solution to the
TSP problem. The algorithm essentially continuously swaps edges in the solution until
the solution no longer improves via this swapping mechanism.
A swap involves removing two edges in the solution, thus leaving two separate paths.
Then two different edges are added to reconnect these two paths. If the new graph has
a better solution than the old one, then it is kept. If not, the change is discarded. Every
possible valid combination of the swapping mechanism is performed.
All of this is just one iteration of the algorithm. We start again with a new solution,
usually better than the previous, and perform every combination of the swapping
mechanism again. We continue iterating this algorithm until either we cannot get an
improved solution, or we can stop it after an arbitrary number of iterations. The runtime
is O(n^3).

## Installing and Running the Program

Compile the program with this command:

```
make tsp
```

Run the program with this command:

```
./tsp [tsp_file_name]
```

Clean files with this command:

```
make clean
```

## Built With

* C++

## Authors

* **Kelvin Lu**
* **Alvin Le**
* **Amanda Grant**
