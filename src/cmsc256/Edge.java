package cmsc256;

public class Edge {
	  private int originVertex;
	  private int destinationVertex;

	  public Edge(int origin, int destination) {
	    this.originVertex = origin;
	    this.destinationVertex = destination;
	  }

	  public int getOriginVertex() {
		return originVertex;
	}

	public void setOriginVertex(int originVertex) {
		this.originVertex = originVertex;
	}

	public int getDestinationVertex() {
		return destinationVertex;
	}

	public void setDestinationVertex(int destinationVertex) {
		this.destinationVertex = destinationVertex;
	}

	public boolean equals(Object o) {
	    return originVertex == ((Edge)o).originVertex && destinationVertex == ((Edge)o).destinationVertex;
	  }
	}
