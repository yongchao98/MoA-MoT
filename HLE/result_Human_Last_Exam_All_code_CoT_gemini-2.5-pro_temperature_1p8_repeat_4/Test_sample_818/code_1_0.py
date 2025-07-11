import itertools

def count_slices(graph_name, num_vertices, edges):
    """
    Counts the number of slices for a given graph G.
    A slice is a partition of edges into two classes where each vertex is incident to at least one edge in each class.
    This is equivalent to finding the number of 2-edge-colorings with no monochromatic vertices, C(G),
    and then N(G) = C(G) / 2.
    """
    num_edges = len(edges)

    # Build an adjacency list where each entry for a vertex is a list of its incident edge indices.
    v_edges = [[] for _ in range(num_vertices)]
    for i, edge in enumerate(edges):
        u, v = edge
        if u < num_vertices and v < num_vertices:
            v_edges[u].append(i)
            v_edges[v].append(i)
        else:
            print(f"Error: vertex index out of range for graph {graph_name}")
            return -1

    valid_colorings = 0
    # Iterate through all 2^num_edges possible colorings.
    # We represent colors with 0 and 1.
    for i in range(2**num_edges):
        coloring = [(i >> j) & 1 for j in range(num_edges)]
        
        is_valid = True
        # Check the condition at each vertex.
        for v_idx in range(num_vertices):
            # Get the colors of edges incident to the current vertex.
            incident_edge_indices = v_edges[v_idx]
            colors = [coloring[idx] for idx in incident_edge_indices]
            
            # Check if all colors are the same (monochromatic vertex).
            if len(set(colors)) == 1:
                is_valid = False
                break
        
        if is_valid:
            valid_colorings += 1
            
    num_slices = valid_colorings // 2
    
    print(f"For {graph_name} (m={num_vertices}): C(G) = {valid_colorings}, N(G) = {num_slices}")
    return num_slices

def solve():
    """
    Determines M(0), M(3), M(5) by checking cubic graphs in increasing order of size.
    """
    print("Starting calculations for M(n)...")

    # M(0): Theoretical result.
    m0 = "none"
    print("\n--- Determining M(0) ---")
    print("For any cubic graph G, the number of slices N(G) is always at least 1.")
    print("Therefore, N(G) cannot be 0, and M(0) does not exist.")

    # M(3):
    m3 = "none"
    print("\n--- Determining M(3) ---")
    # m = 4: The complete graph K4
    k4_v = 4
    k4_e = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    n_k4 = count_slices("K4", k4_v, k4_e)
    if n_k4 % 3 == 0:
        m3 = k4_v
        print(f"N(K4) is a multiple of 3. Since m=4 is the smallest size for a cubic graph, M(3) = {m3}.")
    
    # M(5):
    m5 = "none"
    print("\n--- Determining M(5) ---")
    # m = 4 was already checked for M(3)
    if n_k4 % 5 == 0:
        m5 = k4_v
        print(f"N(K4) is a multiple of 5. Since m=4 is the smallest size, M(5) = {m5}.")
    else:
        print("N(K4) is not a multiple of 5. Checking m=6 graphs...")
        # m = 6: Prism graph and K_3,3
        prism_v = 6
        prism_e = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
        n_prism = count_slices("Prism(6)", prism_v, prism_e)
        
        k33_v = 6
        k33_e = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]
        n_k33 = count_slices("K_3,3", k33_v, k33_e)
        
        if n_prism % 5 == 0 or n_k33 % 5 == 0:
            m5 = 6
            print(f"Found a graph with m=6 whose N(G) is a multiple of 5. M(5) = {m5}.")
        else:
            print("No m=6 graph found. Checking m=8 graphs...")
            # m = 8: Cube graph Q3
            q3_v = 8
            q3_e = []
            for i in range(q3_v):
                for j in range(i + 1, q3_v):
                    # Edge if Hamming distance is 1
                    if bin(i ^ j).count('1') == 1:
                        q3_e.append((i, j))
            n_q3 = count_slices("Cube Q3", q3_v, q3_e)
            
            if n_q3 % 5 == 0:
                m5 = 8
                print(f"N(Q3) is a multiple of 5. M(5) = {m5}.")

    print("\n--- Summary ---")
    print(f"M(0) = {m0}")
    print(f"M(3) = {m3}")
    print(f"M(5) = {m5}")
    print(f"<<<{m0},{m3},{m5}>>>")

if __name__ == '__main__':
    solve()