/********************
** Author: Kelvin Lu, Alvin Le, Amanda Grant
** Date: 03/04/2018
** Description: Main file for final project
********************/

#include <climits>
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <stack>

/********************
** Name: City structure
** Description: Holds the city's id, x-coordinate, y-coordinate, and edge list
********************/
struct city
{
    int id;
    int x;
    int y;
};


/********************
** Name: Distance calculator
** Description: Determine the distance between two cities
********************/
int distance(city* c1, city* c2)
{
    return round(sqrt(pow(c1->x - c2->x, 2) + pow(c1->y - c2->y, 2)));
}

/********************
** Name: Distance matrix creator
** Description: Creates a distance matrix from a list of cities
********************/
std::vector<std::vector<int>> distMatrix(const std::vector<city*>& cities)
{
    int size = cities.size();
    std::vector<std::vector<int>> matrix(size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrix[i].push_back(distance(cities[i], cities[j]));
        }
    }

    return matrix;
}

/********************
** Name: MST-PRIM
** Description: Prim's algorithm for building a minimum spanning tree
** Code Reference: https://www.geeksforgeeks.org/greedy-algorithms-set-5-prims-minimum-spanning-tree-mst-2/
********************/
// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int minKey(std::vector<int> key, std::vector<bool> mstSet, int size)
{
    // Initialize min value
    int min = INT_MAX, min_index;

    for (int v = 0; v < size; v++)
        if (mstSet[v] == false && key[v] < min)
            min = key[v], min_index = v;

    return min_index;
}

// Function to construct and print MST for a graph represented using adjacency
// matrix representation
std::vector<std::vector<int>> primMST(const std::vector<std::vector<int>>& distMatrix)
{
    int size = distMatrix.size();
    std::vector<int> parent (size, -1); // Array to store constructed MST
    std::vector<int> key (size, INT_MAX);   // Key values used to pick minimum weight edge in cut
    std::vector<bool> mstSet (size, false);  // To represent set of vertices not yet included in MST

    // Always include first 1st vertex in MST.
    key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
    parent[0] = -1; // First node is always root of MST

    for (int i = 0; i < size - 1; i++)
    {
        // Pick the minimum key vertex from the set of vertices
        // not yet included in MST
        int u = minKey(key, mstSet, size);

        // Add the picked vertex to the MST Set
        mstSet[u] = true;

        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < size; v++)

            // mstSet[v] is false for vertices not yet included in MST
            // Update the key only if graph[u][v] is smaller than key[v]
            if (distMatrix[u][v] && mstSet[v] == false && distMatrix[u][v] <  key[v])
                parent[v] = u, key[v] = distMatrix[u][v];
    }

    //Add edges from parent into MST adjacency list:
    std::vector<std::vector<int>>mst(parent.size());
    for (int i = 1; i < parent.size(); i++)
    {
        //If node has a parent, add an edge to both vertices in the list:
        mst[i].push_back(parent[i]);
        mst[parent[i]].push_back(i);
    }

    //Return the adjacency list
    return mst;
}

/******************/


/********************
** Name: Matching algorithm
** Description: Algorithm for matching odd degree vertices in the MST to the ones closest to it. Will provide a near-minimal perfect matching
********************/
void match(std::vector<std::vector<int>>& mst, const std::vector<std::vector<int>>& distMatrix)
{
    int mstSize = mst.size();
    //Create a new vector with the indexes of verticees with odd degrees
    std::vector<int> oddVertices;
    for (int i = 0; i < mstSize; i++)
    {
        if (mst[i].size() % 2 != 0)
        {
            oddVertices.push_back(i);
        }
    }

    //Match each odd degree vertex with the odd degree vertex closest to it
    while(!oddVertices.empty())
    {
        int size = oddVertices.size();
        int dist = INT_MAX;
        int close;

        //Find the closest vertex to the current first vertex in the oddVertices vector
        for (int i = 1; i < size; i++)
        {
            if (distMatrix[oddVertices[0]][oddVertices[i]] < dist)
            {
                dist = distMatrix[oddVertices[0]][oddVertices[i]];
                close = i;
            }
        }

        //After match is found, add an edge between them
        mst[oddVertices[0]].push_back(oddVertices[close]);
        mst[oddVertices[close]].push_back(oddVertices[0]);

        //Remove the matched indices
        oddVertices.erase(oddVertices.begin() + close);
        oddVertices.erase(oddVertices.begin());
    }
}


/********************
** Name: Euler Circuit algorithm
** Description: Given an adjacency list, computes the euler circuit for the list using Heirholzer's Algorithm
** Code Reference: https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/ modified for an undirected graph
********************/
std::vector<int> eulerTour(std::vector<std::vector<int> > adj)
{
    // adj represents the adjacency list of
    // the directed graph
    // edge_count represents the number of edges
    // emerging from a vertex
    std::unordered_map<int, int> edge_count;

    for (int i = 0; i<adj.size(); i++)
    {
        //find the count of edges to keep track
        //of unused edges
        edge_count[i] = adj[i].size();
    }

    // Maintain a stack to keep vertices
    std::stack<int> curr_path;

    // vector to store final circuit
    std::vector<int> circuit;

    // start from any vertex
    curr_path.push(0);
    int curr_v = 0; // Current vertex

    while (!curr_path.empty())
    {
        // If there's remaining edge
        if (edge_count[curr_v])
        {
            // Push the vertex
            curr_path.push(curr_v);

            // Find the next vertex using an edge
            int next_v = adj[curr_v].back();

            // and remove that edge from both vertices edge list
            edge_count[curr_v]--;
            adj[curr_v].pop_back();
            for (int i = 0; i < adj[next_v].size(); i++)
            {
                if (adj[next_v][i] == curr_v)
                {
                    adj[next_v].erase(adj[next_v].begin() + i);
                    edge_count[next_v]--;
                    break;
                }
            }

            // Move to next vertex
            curr_v = next_v;
        }

        // back-track to find remaining circuit
        else
        {
            circuit.push_back(curr_v);

            // Back-tracking
            curr_v = curr_path.top();
            curr_path.pop();
        }
    }

    return circuit;
}

/********************
** Name: TSP tour from Euler tour algorithm
** Description: Given an euler circuit, converts it to a hamiltonian tour by deleting redundancies and converting it to a Hamiltonian circuit (Possible due to triangle inequality)
********************/
std::vector<int> hamiltonianTour(std::vector<int> eulerTour)
{
    //Create a visited array corresponding to the vertices
    int size = eulerTour.size();
    std::vector<bool> visited(size , false);
    std::vector<int> tspCircuit;

    for (int i = 0; i< size; i++)
    {
        //If not yet visited, add to tspCircuit
        if (!visited[eulerTour[i]])
        {
            tspCircuit.push_back(eulerTour[i]);
        }
        //Set the vertex to visited
        visited[eulerTour[i]] = true;
    }

    return tspCircuit;
}

/********************
 ** Name: Calculate tour distance. 2-opt helper function.
 ** Description: Calculate total distance of TSP solution.
 *******************/
int calculateTourDistance(std::vector<city*> tour){
    int totalDistance = 0;
    int size = tour.size();

    for(int i = 0; i < size - 1; i++)
    {
        totalDistance += distance(tour[i], tour[i+1]);
    }

    return totalDistance;
}

/********************
** Name: Driver Function
** Description: Input a vector of cities and return the TSP tour in a vector
********************/
std::vector<city*> driver(const std::vector<city*> &cities, int &cost)
{
    //Create distnace matrix from cities
    std::vector<std::vector<int>> matrix = distMatrix(cities);

    //Create minimum spanning tree adjacency list
    std::vector<std::vector<int>> adjacencyList = primMST(matrix);

    //Match odd vertices in the minimum spanning tree and add those edges to the MST
    match(adjacencyList, matrix);

    //Find a euler cycle in the resulting adjacency list
    std::vector<int> euler = eulerTour(adjacencyList);

    //Convert euler cycle to hamiltonian cycle
    std::vector<int> hamiltonian = hamiltonianTour(euler);

    //Build the TSP tour from the hamiltonian cycle
    std::vector<city*> tspTour;
    for (int i = 0; i < hamiltonian.size() - 1; i++)
    {
        cost += matrix[hamiltonian[i]][hamiltonian[i + 1]];
        tspTour.push_back(cities[hamiltonian[i]]);
    }

    //Add last vertex and edge to end cycle
    int last = hamiltonian.back();
    tspTour.push_back(cities[last]);
    cost += matrix[hamiltonian[0]][last];

    return tspTour;
}



/******************
 ** Name: 2-opt swap
 ** Description: Helper function for TwoOpt. Performs edge swap in tour.
 ** Returns: New tour that resulted from swapping.
 ******************/
std::vector<city*> TwoOptSwap(std::vector<city*> tour, int i, int k){
    std::vector<city*> newTour;
    int size = tour.size();

    //1. add tour[0] to tour[i-1] in order to new tour.
    for(int j = 0; j <= i - 1; j++)
    {
        newTour.push_back(tour[j]);
    }

    //2. add tour[i] to tour[k] in reverse order to new tour.
    for (int j = k; j >= i; j--)
    {
        newTour.push_back(tour[j]);
    }

    //3. add tour[k+1] to end in order to new tour.
    for (int j = k + 1; j < size; j++)
    {
        newTour.push_back(tour[j]);
    }
    return newTour;
}

/*******************
 ** Name: TwoOpt optimization algorithm.
 ** Description: Optimizes a solution to the TSP problem.
 **    Swaps edges in the solution to make improvements until no further improvements can be made.
 ** Recieves: A solution to the TSP problem.
 ** Returns: Returns the total distance of the new solution. Also changes the input tour to match the new solution.
 ******************/
int TwoOpt(std::vector<city*> &tour){

    //Add first city to the end to complete the cycle.
    tour.push_back(tour[0]);

    int numCities = tour.size() - 1;
    int solution = calculateTourDistance(tour);
    bool improve;

    //This do-while loops swaps edges in the solution until no improvement can be made via the swapping mechanism.
    do
    {
        improve = false;
        for (int i = 1; i < numCities - 1; i++)
        {
            for (int k = i + 1; k < numCities; k++)
            {
                //Perform swap and obtain a new tour and solution.
                std::vector<city*> newTour = TwoOptSwap(tour, i, k);
                int newSolution = calculateTourDistance(newTour);


                //Check if the new solution is better than the old solution.
                //If yes, update the tour and solution.
                //Else, keep old tour and solution.
                if (newSolution < solution){\
                    improve = true;
                    tour = newTour;
                    solution = newSolution;
                }
            }
        }
    }while(improve == true);

    //Remove last city to obtain proper solution format.
    tour.pop_back();

    return solution;
}


/********************
** Name: File I/O
** Description: Accepts filename through command line, parses file for city information, and outputs results to a new file
********************/
int main(int argc, char **argv)
{
    std::ifstream ifs;
    std::ofstream ofs;
    std::string curr;
    clock_t time1, time2;
    char * cstr;
    char * token;
    int cost = 0;

    std::vector<city*> cities;

    std::string filename = argv[1];
    ifs.open(filename);

    //Check if datafile exists
    if (!ifs)
    {
        std::cout << "Error: File not found." << std::endl;
        ifs.close();
        return 1;
    }

    //While end of file has not been reached
    while(getline(ifs, curr))
    {
        //Create a new city
        city* c = new city;

        //Read and tokenize the string
        cstr = new char[curr.length() + 1];
        strcpy(cstr, curr.c_str());
        token = strtok(cstr, " \r\n");

        //Set city values to inputs
        c->id = std::stoi(token);
        c->x = std::stoi(strtok(NULL, " \r\n"));
        c->y = std::stoi(strtok(NULL, " \r\n"));

        //Add city to cities vector
        cities.push_back(c);

        //Deallocate c-string
        delete[] cstr;

    }

    //Close input file stream
    ifs.close();

    //Get start time.
    time1 = clock();

    //Find tour
    std::vector<city*> tour = driver(cities, cost);

    //Run 2-opt if the size of the tour is less than 1000.
    //Otherwise, the run will take too much time.
    if(tour.size() < 1000)
    {
        //Perform 2-opt optimization algorithm.
        cost = TwoOpt(tour);
    }

    //Get finish time.
    time2 = clock();

    //Calculate and print run time.
    float total_time = float(time2) - float(time1);
    float seconds = total_time / CLOCKS_PER_SEC;
    std::cout << "Runtime: " <<  seconds << " seconds" << std::endl;


    //Open output file stream and output tour.
    ofs.open(filename + ".tour");

    ofs << cost << "\n";
    for (int i = 0; i < tour.size(); i++)
    {
        ofs << tour[i]->id << "\n";
    }

    ofs.close();

    //Free allocated memory
    for (int i = 0; i < cities.size(); i++)
    {
        delete cities[i];
    }

}
