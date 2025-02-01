# UnweightedGraph Implementation
This repository contains an implementation of an **UnweightedGraph** for undirected graphs in Java. Models a graph with various functionalities such as adding vertices and edges and traversing the graph using DFS and BFS.
## Features
### Core Functionalities
1. **Vertices and Edges:**
    - Add and remove graph vertices and edges.
    - Store vertices and edges using an adjacency list.

2. **Graph Construction:**
    - Construct graphs using:
        - Arrays or Lists of vertices and edges
        - Integer vertices with edges defined by arrays or lists.

    - Load and build graphs from input files using the `graphFromFile` method.

3. **Graph Traversal:**
    - Perform **Depth-First Search (DFS)** and **Breadth-First Search (BFS)**, returning search trees that include information about search order and parent relationships.

### Graph Utility Methods
- **Graph Properties:**
    - `isComplete`: Check if the graph is complete (i.e., every vertex is connected to all other vertices).
    - `isConnected`: Determine if the graph is connected.
    - `hasCycle`: Check if the graph contains a cycle.

- **Adjacency and Paths:**
    - `areAdjacent`: Check if two vertices are directly connected.
    - `getShortestPath`: Find the shortest path between two vertices using BFS.

- **Degree and Neighbors:**
    - Get the degree of a vertex and list its neighboring vertices.

- **Print Graph Data:**
    - Display graph edges and traversal tree details.

## Program Overview
### `UnweightedGraph<V>`
The primary class that implements the `Graph<V>` interface and provides the following features:
1. **Internal Structure:**
    - Stores vertices (`List<V>`) and edges (`List<List<Edge>>`) using an adjacency list.

2. **Key Methods:**
    - `addVertex(V vertex)`: Add a vertex to the graph.
    - `addEdge(int origin, int destination)`: Add a directed edge between vertices.
    - `getNeighbors(int index)`: Return a list of neighboring vertices for a specific vertex.
    - `getDegree(int vertex)`: Get the degree of a vertex.
    - Graph file parsing via `graphFromFile(String fileName)` to construct a graph from a `.txt file format`.

3. **Traversal Methods:**
    - `getDepthFirstSearchTree(int startingVertex)` and `getBreadthFirstSearchTree(int vertex)` return `SearchTree` objects representing the traversal.

4. **Graph Analysis:**
    - Methods like `isComplete`, `isConnected`, and `hasCycle` (defined above).

### `Edge`
The `Edge` class models a graph edge with origin and destination vertices:
- `Edge(int origin, int destination)`: Constructor to define the start and end vertices of an edge.
- Getter and setter methods for accessing origin and destination.
- Overridden `equals` method for comparison.

### Inner Class: `SearchTree`
The `SearchTree` class represents the result of a DFS or BFS traversal of the graph:
- Stores:
    - Root vertex.
    - Parent of each vertex in the tree.
    - Search order of vertices.

- Methods include:
    - `getPath(int index)`: Retrieve the path from a given vertex to the root.
    - `printPath(int index)`: Display the path from the root to a given vertex.
    - `printTree()`: Print the entire traversal tree.
