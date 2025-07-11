import itertools

def count_slices(graph_name, graph_adj):
    """
    Counts the number of slices for a given cubic graph.

    A slice is a 2-edge-coloring where no vertex is monochromatic.
    The number of slices N(G) is half the number of such colorings.
    """
    vertices = list(graph_adj.keys())
    
    # Create a stable, sorted list of edges to assign indices
    edges = set()
    for u in vertices:
        for v in graph_adj[u]:
            edges.add(tuple(sorted((u, v))))
    
    sorted_edges = sorted(list(edges))
    num_edges = len(sorted_edges)
    edge_map = {edge: i for i, edge in enumerate(sorted_edges)}

    # Map each vertex to the indices of its incident edges
    incident_edge_indices = {v: [] for v in vertices}
    for v in vertices:
        for neighbor in graph_adj[v]:
            edge = tuple(sorted((v, neighbor)))
            incident_edge_indices[v].append(edge_map[edge])

    valid_coloring_count = 0
    # Iterate through all 2^num_edges possible edge colorings
    for i in range(2**num_edges):
        is_valid_coloring = True
        for v in vertices:
            # Get the colors of the 3 edges incident to vertex v
            colors = [
                (i >> edge_idx) & 1 for edge_idx in incident_edge_indices[v]
            ]
            # Check if all colors are the same (monochromatic vertex)
            if colors[0] == colors[1] == colors[2]:
                is_valid_coloring = False
                break
        
        if is_valid_coloring:
            valid_coloring_count += 1
            
    # The number of partitions (slices) is half the number of valid colorings,
    # as swapping colors gives the same partition.
    num_slices = valid_coloring_count // 2
    return num_slices

# --- Graph Definitions ---

# K4: Complete graph on 4 vertices (m=4)
K4 = {
    0: [1, 2, 3],
    1: [0, 2, 3],
    2: [0, 1, 3],
    3: [0, 1, 2]
}

# K(3,3): Complete bipartite graph (Utility graph) (m=6)
K3_3 = {
    0: [3, 4, 5], 1: [3, 4, 5], 2: [3, 4, 5],
    3: [0, 1, 2], 4: [0, 1, 2], 5: [0, 1, 2]
}

# Y3: Prism graph (C3 x K2) (m=6)
Y3 = {
    0: [1, 2, 3],
    1: [0, 2, 4],
    2: [0, 1, 5],
    3: [0, 4, 5],
    4: [1, 3, 5],
    5: [2, 3, 4]
}

# --- Main Logic and Output ---

# Determine M(0)
print("Analysis for M(0):")
print("A 'slice' corresponds to a decomposition of the graph's edges into two spanning subgraphs")
print("where each vertex has at least one incident edge from each subgraph.")
print("It is a known result that every cubic graph has such a slice. Therefore, N(G) >= 1 for any cubic graph G.")
print("For N(G) to be a multiple of 0, N(G) must be 0, which is not possible.")
m0 = "none"
print(f"Result for M(0): {m0}\n")

# Determine M(3)
print("Analysis for M(3):")
print("We need the smallest vertex count 'm' for a cubic graph G such that N(G) is a multiple of 3.")
print("The smallest cubic graph is K4 (m=4). Let's calculate N(K4).")
n_k4 = count_slices("K4", K4)
print(f"The number of slices for K4 is N(K4) = {n_k4}.")
if n_k4 % 3 == 0:
    m3 = 4
    print(f"Since {n_k4} is a multiple of 3, the smallest m is 4.")
else:
    # This branch is not expected to be taken
    m3 = "undetermined" 
    print(f"N(K4) is not a multiple of 3. Need to check larger graphs.")
print(f"Result for M(3): {m3}\n")


# Determine M(5)
print("Analysis for M(5):")
print("We need the smallest 'm' such that N(G) is a multiple of 5.")
print("First, check m=4: N(K4) = 9.")
if n_k4 % 5 == 0:
    m5 = 4
    print(f"Since {n_k4} is a multiple of 5, the smallest m is 4.")
else:
    print(f"N(K4) = {n_k4}, which is not a multiple of 5. We check the next possible size, m=6.")
    print("There are two cubic graphs with 6 vertices: the Prism graph (Y3) and K(3,3).")
    
    n_y3 = count_slices("Y3", Y3)
    print(f"The number of slices for the Prism graph is N(Y3) = {n_y3}.")
    
    n_k3_3 = count_slices("K(3,3)", K3_3)
    print(f"The number of slices for K(3,3) is N(K(3,3)) = {n_k3_3}.")

    if n_y3 % 5 == 0 or n_k3_3 % 5 == 0:
        m5 = 6
        print(f"One of these, N(K(3,3))={n_k3_3}, is a multiple of 5. So the smallest m is 6.")
    else:
        # This branch is not expected to be taken
        m5 = "undetermined"
        print("Neither is a multiple of 5. Need to check larger graphs.")

print(f"Result for M(5): {m5}\n")

print("---")
print("Final Combined Answer: M(0),M(3),M(5)")
print(f"{m0},{m3},{m5}")