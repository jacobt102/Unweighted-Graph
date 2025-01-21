package cmsc256;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class UnweightedGraph<V> implements Graph<V> {
	private List<V> vertices = new ArrayList<>(); // Store vertices
	private List<List<Edge>> neighbors = new ArrayList<>(); // Adjacency lists

	/** Construct an empty graph */
	public UnweightedGraph() {
	}

	/** Construct a graph from vertices and edges stored in arrays */
	public UnweightedGraph(V[] vertices, int[][] edges) {
		for (int i = 0; i < vertices.length; i++)
			addVertex(vertices[i]);

		createAdjacencyLists(edges, vertices.length);
	}

	/** Construct a graph from vertices and edges stored in List */
	public UnweightedGraph(List<V> vertices, List<Edge> edges) {
		for (int i = 0; i < vertices.size(); i++)
			addVertex(vertices.get(i));

		createAdjacencyLists(edges, vertices.size());
	}

	/** Construct a graph for integer vertices 0, 1, 2 and edge list */
	public UnweightedGraph(List<Edge> edges, int numberOfVertices) {
		for (int i = 0; i < numberOfVertices; i++) 
			addVertex((V)(Integer.valueOf(i))); // vertices are {0, 1, ...}

		createAdjacencyLists(edges, numberOfVertices);
	}

	/** Construct a graph from integer vertices 0, 1, and edge array */
	public UnweightedGraph(int[][] edges, int numberOfVertices) {
		for (int i = 0; i < numberOfVertices; i++) 
			addVertex((V)(Integer.valueOf(i))); // vertices are {0, 1, ...}

		createAdjacencyLists(edges, numberOfVertices);
	}

	/** Create adjacency lists for each vertex */
	private void createAdjacencyLists(
			int[][] edges, int numberOfVertices) {
		for (int i = 0; i < edges.length; i++) {
			addEdge(edges[i][0], edges[i][1]);
		}
	}

	/** Create adjacency lists for each vertex */
	private void createAdjacencyLists(
			List<Edge> edges, int numberOfVertices) {
		for (Edge edge: edges) {
			addEdge(edge.getOriginVertex(), edge.getDestinationVertex());
		}
	}

	@Override /** Return the number of vertices in the graph */
	public int getSize() {
		return vertices.size();
	}

	@Override /** Return the vertices in the graph */
	public List<V> getVertices() {
		return vertices;
	}

	@Override /** Return the object for the specified vertex */
	public V getVertex(int index) {
		return vertices.get(index);
	}

	@Override /** Return the index for the specified vertex object */
	public int getIndex(V v) {
		return vertices.indexOf(v);
	}

	@Override /** Return the neighbors of the specified vertex */
	public List<Integer> getNeighbors(int index) {
		List<Integer> result = new ArrayList<>();
		for (Edge e: neighbors.get(index))
			result.add(e.getDestinationVertex());

		return result;
	}

	@Override /** Return the degree for a specified vertex */
	public int getDegree(int v) {
		return neighbors.get(v).size();
	}

	@Override /** Print the edges */
	public void printEdges() {
		for (int u = 0; u < neighbors.size(); u++) {
			System.out.print(getVertex(u) + " (" + u + "): ");
			for (Edge e: neighbors.get(u)) {
				System.out.print("(" + getVertex(e.getOriginVertex()) + ", " +
						getVertex(e.getDestinationVertex()) + ") ");
			}
			System.out.println();
		}
	}

	@Override /** Clear the graph */
	public void clear() {
		vertices.clear();
		neighbors.clear();
	}

	@Override /** Add a vertex to the graph */  
	public boolean addVertex(V vertex) {
		if (!vertices.contains(vertex)) {
			vertices.add(vertex);
			neighbors.add(new ArrayList<Edge>());
			return true;
		}
		else {
			return false;
		}
	}

	@Override /** Add an edge to the graph */  
	public boolean addEdge(Edge e) {
		if (e.getOriginVertex() < 0 || e.getOriginVertex() > getSize() - 1)
			throw new IllegalArgumentException("No such index: " + e.getOriginVertex());

		if (e.getDestinationVertex() < 0 || e.getDestinationVertex() > getSize() - 1)
			throw new IllegalArgumentException("No such index: " + e.getDestinationVertex());

		if (!neighbors.get(e.getOriginVertex()).contains(e)) {
			neighbors.get(e.getOriginVertex()).add(e);
			return true;
		}
		else {
			return false;
		}
	}

	@Override /** Add an edge to the graph */  
	public boolean addEdge(int u, int v) {
		return addEdge(new Edge(u, v));
	}

@Override /** Obtain a DFS tree starting from vertex v */
	public SearchTree getDepthFirstSearchTree(int startingVertex) {
		List<Integer> searchOrder = new ArrayList<>();
		int[] parent = new int[vertices.size()];
		for (int i = 0; i < parent.length; i++)
			parent[i] = -1; // Initialize parent[i] to -1

		// Mark visited vertices
		boolean[] isVisited = new boolean[vertices.size()];

		// Recursively search
		dfs(startingVertex, parent, searchOrder, isVisited);

		// Return a search tree
		return new SearchTree(startingVertex, parent, searchOrder);
	}

	/** Recursive method for DFS search */
	private void dfs(int v, int[] parent, List<Integer> searchOrder, boolean[] isVisited) {
		// Store the visited vertex
		searchOrder.add(v);
		isVisited[v] = true; // Vertex v visited

		for (Edge e : neighbors.get(v)) { // Note that e.u is v
			if (!isVisited[e.getDestinationVertex()]) { // e.v is w in Listing 28.8
				parent[e.getDestinationVertex()] = v; // The parent of w is v
				dfs(e.getDestinationVertex(), parent, searchOrder, isVisited); // Recursive search
			}
		}
	}

	@Override /** Starting bfs search from vertex v */
	public SearchTree getBreadthFirstSearchTree(int vertex) {
		List<Integer> searchOrder = new ArrayList<>();
		int[] parent = new int[vertices.size()];
		for (int i = 0; i < parent.length; i++)
			parent[i] = -1; // Initialize parent[i] to -1

		// list used as a queue
		java.util.LinkedList<Integer> queue = new java.util.LinkedList<>(); 
		boolean[] isVisited = new boolean[vertices.size()];
		queue.offer(vertex); 		// Enqueue vertex
		isVisited[vertex] = true; 	// Mark it visited

		while (!queue.isEmpty()) {
			int u = queue.poll(); 	// Dequeue to u
			searchOrder.add(u); 		// u searched
			for (Edge e: neighbors.get(u)) { // Note that e.u is u
				if (!isVisited[e.getDestinationVertex()]) { 
					queue.offer(e.getDestinationVertex()); 		// Enqueue w
					parent[e.getDestinationVertex()] = u; 		// The parent of w is u
					isVisited[e.getDestinationVertex()] = true; 	// Mark w visited
				}
			}
		}

		return new SearchTree(vertex, parent, searchOrder);
	}
	
	
	/** SearchTree inner class inside the UnweightedGraph class */
	public class SearchTree {
		private int root; // The root of the tree
		private int[] parent; // Store the parent of each vertex
		private List<Integer> searchOrder; // Store the search order

		/** Construct a tree with root, parent, and searchOrder */
		public SearchTree(int root, int[] parent, 
				List<Integer> searchOrder) {
			this.root = root;
			this.parent = parent;
			this.searchOrder = searchOrder;
		}

		/** Return the root of the tree */
		public int getRoot() {
			return root;
		}

		/** Return the parent of vertex v */
		public int getParent(int v) {
			return parent[v];
		}

		/** Return an array representing search order */
		public List<Integer> getSearchOrder() {
			return searchOrder;
		}

		/** Return number of vertices found */
		public int getNumberOfVerticesFound() {
			return searchOrder.size();
		}

		/** Return the path of vertices from a vertex to the root */
		public List<V> getPath(int index) {
			ArrayList<V> path = new ArrayList<>();

			do {
				path.add(vertices.get(index));
				index = parent[index];
			}
			while (index != -1);

			return path;
		}

		/** Print a path from the root to vertex v */
		public void printPath(int index) {
			List<V> path = getPath(index);
			System.out.print("A path from " + vertices.get(root) + " to " +
					vertices.get(index) + ": ");
			for (int i = path.size() - 1; i >= 0; i--)
				System.out.print(path.get(i) + " ");
		}

		/** Print the whole tree */
		public void printTree() {
			System.out.println("Root is: " + vertices.get(root));
			System.out.print("Edges: ");
			for (int i = 0; i < parent.length; i++) {
				if (parent[i] != -1) {
					// Display an edge
					System.out.print("(" + vertices.get(parent[i]) + ", " +
							vertices.get(i) + ") ");
				}
			}
			System.out.println();
		}
	}
	/*  End implementation of inner SearchTree class */

        /*  Lab methods to be implemented	*/
	
	/**
	 * Reads the input file to extract the graph data
	 * @param fileName      The String name of the file
	 * @return              UnweightedGraph of data from the File
	 * @throws FileNotFoundException, NumberFormatException
	 */
	public static UnweightedGraph<Integer> graphFromFile(String fileName) throws FileNotFoundException{
		File file = new File(fileName);
		if(!file.canRead()){
			throw new FileNotFoundException();
		}
		Scanner readFile = new Scanner(file);
		int numVertices = Integer.parseInt(readFile.nextLine());
		UnweightedGraph<Integer> returnGraph = new UnweightedGraph<>();
		for(int i = 0; i < numVertices; i++){
			returnGraph.addVertex(i);
		}
		while(readFile.hasNextLine()){
			String[] nums = readFile.nextLine().split(" ");
			for (int i = 1; i < nums.length; i++) {
				returnGraph.addEdge(Integer.parseInt(nums[0]),Integer.parseInt(nums[i]));
			}
		}

		return returnGraph;
	}

	/**
	 * @return 		true if this graph is complete
	 */
	public boolean isComplete() {
        return getDegree(1) == getSize() - 1;
	}
	
	/**
	 * @param origin 	the origin vertex
	 * @param destination	the destination vertex
	 * @return 		true if the destination vertex is adjacent 
	 * 			to the origin vertex; false otherwise
	 */
	public boolean areAdjacent(int origin, int destination) {
		List<Integer> list = getNeighbors(destination);
        return list.contains(origin);
	}
	
	/**
	 * @return 	true if this graph is connected
	 */
	public boolean isConnected() {
		printEdges();
		SearchTree tree = getDepthFirstSearchTree(0);
		if(getSize() == 1){
			return false;
		}
        return tree.getNumberOfVerticesFound() == getSize();
	}



	/**
	 * @return 	true if this graph has a cycle
	 */
	public boolean hasCycle() {
		return false;
	}
	/**
	 * @param origin	the origin vertex
	 * @param destination	the destination vertex
	 * @return		A List containing the shortest path from 
	 * 			the origin vertex to the destination vertex
	 */
	public List<V> getShortestPath(int origin, int destination){
		if(vertices.isEmpty() || getVertex(destination)==null || getVertex(origin) == null){
			return null;
		}
		SearchTree tree = getBreadthFirstSearchTree(destination);

		List<V> returnList = tree.getPath(origin);
		if(returnList.size()>1){
			return returnList;
		}
		else{
			return null;
		}
	}

	public static void main(String[] args) {
		try {
			graphFromFile("test1.txt").getShortestPath(0,1);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
	}
	


}

