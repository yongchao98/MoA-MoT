import math

def count_2_matchings(graph_adj):
    """
    Calculates the number of 2-matchings in a graph.

    A 2-matching is a set of 2 edges that do not share any vertices.
    This function implements a formula that runs in O(|V| + |E|) time,
    which is subcubic in the number of vertices |V|.

    The formula is:
    Num_2_matchings = (Total Pairs of Edges) - (Pairs of Incident Edges)
                    = C(|E|, 2) - sum(C(deg(v), 2) for v in V)

    Args:
        graph_adj: The graph represented as an adjacency list (dictionary).
    """
    num_vertices = len(graph_adj)
    if num_vertices == 0:
        print("Graph is empty.")
        return

    degrees = [len(graph_adj[v]) for v in range(num_vertices)]
    
    # Calculate total number of edges, |E|
    # sum(degrees) is 2*|E|
    num_edges = sum(degrees) // 2
    
    # Calculate total number of edge pairs: C(|E|, 2)
    # Using C(n,k) = n*(n-1)//2 for k=2
    if num_edges < 2:
        total_edge_pairs = 0
    else:
        total_edge_pairs = num_edges * (num_edges - 1) // 2
        
    # Calculate number of incident edge pairs
    # This is sum over all vertices v of C(degree(v), 2)
    incident_edge_pairs = 0
    for deg in degrees:
        if deg >= 2:
            incident_edge_pairs += deg * (deg - 1) // 2
            
    # The final count for 2-matchings
    num_2_matchings = total_edge_pairs - incident_edge_pairs

    print("--- Calculating the number of 2-matchings ---")
    print("The final equation is: Num_2_matchings = Total_Edge_Pairs - Incident_Edge_Pairs")
    print("\n--- Intermediate values for the equation ---")
    print(f"Number of edges |E|: {num_edges}")
    print(f"Total pairs of edges C(|E|, 2): {total_edge_pairs}")
    print(f"Number of incident edge pairs sum(C(deg(v), 2)): {incident_edge_pairs}")
    
    print("\n--- Final Result ---")
    print(f"The number of 2-matchings is: {num_2_matchings}")


if __name__ == '__main__':
    # Example: A pentagon graph C5
    # It has 5 vertices and 5 edges. Each vertex has degree 2.
    # The number of 2-matchings should be 5.
    c5_graph = {
        0: [1, 4],
        1: [0, 2],
        2: [1, 3],
        3: [2, 4],
        4: [0, 3]
    }
    count_2_matchings(c5_graph)
