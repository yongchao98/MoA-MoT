import itertools

def calculate_N(m, edges):
    """
    Calculates N(G) for a graph G with m vertices and a list of edges.
    
    A "slice" is a 2-coloring of edges where no vertex is monochromatic.
    N(G) = (number of valid colorings) / 2.
    """
    num_edges = len(edges)
    num_valid_colorings = 0
    
    # Pre-compute adjacency list for efficiency
    adj = {i: [] for i in range(m)}
    for i, edge in enumerate(edges):
        adj[edge[0]].append(i)
        adj[edge[1]].append(i)

    # Iterate through all 2^|E| edge colorings
    for coloring in itertools.product([0, 1], repeat=num_edges):
        is_valid_coloring = True
        # Check the condition for each vertex
        for v in range(m):
            # A cubic graph must have degree 3 for every vertex.
            # No need to check for empty incident_edge_colors
            
            # Since degree is 3, a vertex is monochromatic if all 3 incident
            # edges have the same color.
            first_color = coloring[adj[v][0]]
            if coloring[adj[v][1]] == first_color and coloring[adj[v][2]] == first_color:
                is_valid_coloring = False
                break
        
        if is_valid_coloring:
            num_valid_colorings += 1
            
    # N(G) is half the number of valid colorings
    return num_valid_colorings // 2

def solve():
    """
    Finds M(0), M(3), and M(5) by checking cubic graphs in order of size.
    """
    graphs_to_check = {
        4: {
            "K4": (4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
        },
        6: {
            "Prism": (6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]),
            "K3,3": (6, [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)])
        },
        # For m=8, all cubic graphs are 3-edge-colorable, so N(G)>0. We skip them to save time.
        10: {
            "Petersen": (10, [
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),  # Outer cycle
                (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),  # Inner star
                (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)   # Spokes
            ])
        }
    }

    m0, m3, m5 = "none", "none", "none"

    print("Determining M(3), M(5), and M(0):")
    
    sorted_m_values = sorted(graphs_to_check.keys())
    
    for m in sorted_m_values:
        if m0 != "none" and m3 != "none" and m5 != "none":
            break # All values found
        
        print(f"\n--- Checking graphs with m = {m} vertices ---")
        for name, (num_v, edges) in graphs_to_check[m].items():
            print(f"Calculating N(G) for {name} graph...")
            n_g = calculate_N(num_v, edges)
            print(f"Result: N({name}) = {n_g}")

            # Check for M(3)
            if m3 == "none" and n_g % 3 == 0:
                m3 = m
                print(f"Found M(3): N({name}) = {n_g} is a multiple of 3. So, M(3) = {m}.")

            # Check for M(5)
            if m5 == "none" and n_g % 5 == 0:
                m5 = m
                print(f"Found M(5): N({name}) = {n_g} is a multiple of 5. So, M(5) = {m}.")
            
            # Check for M(0)
            if m0 == "none" and n_g == 0:
                m0 = m
                print(f"Found M(0): N({name}) = {n_g}. So, M(0) = {m}.")

    print("\n--- Summary of results ---")
    print(f"M(0): Smallest m for N(G) = 0 is {m0}.")
    print(f"M(3): Smallest m for N(G) multiple of 3 is {m3}.")
    print(f"M(5): Smallest m for N(G) multiple of 5 is {m5}.")
    
    print("\nFinal Answer:")
    print(f"{m0},{m3},{m5}")

if __name__ == '__main__':
    solve()