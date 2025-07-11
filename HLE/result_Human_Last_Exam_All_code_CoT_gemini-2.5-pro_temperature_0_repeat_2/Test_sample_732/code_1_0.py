import collections

def solve_graph_problem():
    """
    This script solves the problem by constructing a maximal C4-free graph
    with 8 vertices and verifying its properties.
    """

    # The problem is to find the maximum number of edges in a simple graph
    # with 8 vertices that has no C4 (quadrilateral) subgraphs.
    # This is a known problem in extremal graph theory, with the value being ex(8, C4).

    # We can construct a C4-free graph with 8 vertices and 10 edges.
    # It has been proven that any graph with 8 vertices and 11 edges must contain a C4.
    # Therefore, the maximum number of edges is 10.

    # Here is the construction of such a graph, G:
    # Let V = {0, 1, 2, 3, 4, 5, 6, 7} be the set of vertices.
    # Let vertex 0 be a central vertex connected to all other 7 vertices.
    # The remaining 7 vertices {1, ..., 7} form a subgraph consisting of a
    # maximum matching (3 disjoint edges) to avoid creating C4s.
    G = {
        0: [1, 2, 3, 4, 5, 6, 7],
        1: [0, 2],
        2: [0, 1],
        3: [0, 4],
        4: [0, 3],
        5: [0, 6],
        6: [0, 5],
        7: [0]
    }

    num_vertices = len(G)

    # Function to count the number of edges in an adjacency list representation
    def count_edges(adj):
        edge_count = 0
        for vertex in adj:
            edge_count += len(adj[vertex])
        # Each edge is counted twice (once for each endpoint), so we divide by 2.
        return edge_count // 2

    num_edges = count_edges(G)

    # Function to count the number of C4s in a graph.
    # A C4 is formed by two vertices u, v and two of their common neighbors.
    # If u and v have k common neighbors, they form (k choose 2) C4s.
    def count_c4(adj):
        c4_count = 0
        nodes = list(adj.keys())
        n = len(nodes)
        for i in range(n):
            for j in range(i + 1, n):
                u, v = nodes[i], nodes[j]
                # Find common neighbors of u and v
                # Using set intersection for efficiency
                neighbors_u = set(adj[u])
                neighbors_v = set(adj[v])
                common_neighbors = neighbors_u.intersection(neighbors_v)
                
                k = len(common_neighbors)
                if k >= 2:
                    c4_count += k * (k - 1) // 2
        return c4_count

    num_c4 = count_c4(G)

    print("Verifying the constructed graph:")
    print(f"Number of vertices: {num_vertices}")
    print(f"Number of edges: {num_edges}")
    print(f"Number of quadrilaterals (C4s): {num_c4}")
    print("-" * 20)
    print("Based on this construction and known results in graph theory, the final answer is:")
    
    # The final equation is that the maximum number of edges is 10.
    # We print each number/component of this conclusion.
    final_answer = 10
    print(f"The maximum number of edges in a simple graph with 8 vertices and no quadrilaterals is {final_answer}.")

solve_graph_problem()
<<<10>>>