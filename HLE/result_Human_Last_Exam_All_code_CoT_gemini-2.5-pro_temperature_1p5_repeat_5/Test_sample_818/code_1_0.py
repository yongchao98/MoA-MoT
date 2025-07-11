import itertools

def count_slices(num_vertices, edges):
    """
    Counts the number of slices for a given graph G=(V,E).
    A slice is a partition of E into two classes such that each vertex is
    incident to at least one edge in each class. This is equivalent to
    finding the number of 2-edge-colorings where no vertex is monochromatic.
    The number of slices N(G) is half the number of such colorings, C(G).
    """
    num_edges = len(edges)
    
    # Create adjacency list where adj[v] = list of edge indices incident to v
    adj = [[] for _ in range(num_vertices)]
    for i, (u, v) in enumerate(edges):
        adj[u].append(i)
        adj[v].append(i)

    valid_colorings = 0
    # Iterate through all 2^|E| possible 2-edge-colorings
    # A coloring is represented by an integer from 0 to 2^|E| - 1
    for i in range(2**num_edges):
        is_valid = True
        # Check each vertex for the slice condition
        for v in range(num_vertices):
            colors = []
            for edge_idx in adj[v]:
                # Get the color of the edge (0 or 1)
                color = (i >> edge_idx) & 1
                colors.append(color)
            
            # A vertex is monochromatic if the sum of incident edge colors is 0 or 3
            if sum(colors) == 0 or sum(colors) == 3:
                is_valid = False
                break
        
        if is_valid:
            valid_colorings += 1
            
    # The number of slices is half the number of valid colorings
    num_slices = valid_colorings / 2
    return int(num_slices)

def solve():
    """
    Finds M(0), M(3), and M(5) by checking cubic graphs in increasing order of size.
    """
    print("Step-by-step calculation:")

    # M(0): Based on the reasoning that no cubic graph G has N(G)=0.
    m_0 = "none"
    print("M(0): Assuming Bouchet's conjecture is true, no cubic graph has N(G)=0. So M(0) is none.")

    # M(3) and M(5)
    m_vals = {3: None, 5: None}
    
    # List of non-isomorphic cubic graphs to check, by (num_vertices, name, edge_list)
    # Edge lists are 0-indexed.
    graphs_to_check = [
        # 4-vertex graphs
        (4, "K4", [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]),
        # 6-vertex graphs
        (6, "Prism(Y3)", [(0,1),(1,2),(2,0), (3,4),(4,5),(5,3), (0,3),(1,4),(2,5)]),
        (6, "K_3,3", [(0,3),(0,4),(0,5), (1,3),(1,4),(1,5), (2,3),(2,4),(2,5)]),
        # 8-vertex graphs (all 5 non-isomorphic ones)
        (8, "Cube(Q3)", [(0,1),(0,2),(0,4),(1,3),(1,5),(2,3),(2,6),(3,7),(4,5),(4,6),(5,7),(6,7)]),
        (8, "G(8,2)", [(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(2,7),(3,4),(3,6),(4,7),(5,6),(5,7)]),
        (8, "G(8,3)", [(0,1),(0,4),(0,5),(1,2),(1,6),(2,3),(2,7),(3,4),(3,7),(4,5),(5,6),(6,7)]),
        (8, "G(8,4)", [(0,1),(0,3),(0,4),(1,2),(1,5),(2,3),(2,6),(3,7),(4,5),(4,7),(5,6),(6,7)]),
        (8, "G(8,5)", [(0,1),(0,3),(0,4),(1,2),(1,5),(2,3),(2,6),(3,7),(4,6),(4,7),(5,6),(5,7)]),
    ]

    for v_count, name, edges in graphs_to_check:
        if m_vals[3] is not None and m_vals[5] is not None:
            break

        n_g = count_slices(v_count, edges)
        print(f"Graph: {name} ({v_count} vertices), N(G) = {n_g}")

        if n_g % 3 == 0 and m_vals[3] is None:
            m_vals[3] = v_count
            print(f"Found M(3)! N({name})={n_g} is a multiple of 3. M(3) = {v_count}")
        
        if n_g % 5 == 0 and m_vals[5] is None:
            m_vals[5] = v_count
            print(f"Found M(5)! N({name})={n_g} is a multiple of 5. M(5) = {v_count}")

    print("\nFinal Answer Calculation Complete.")
    print(f"M(0) = {m_0}")
    print(f"M(3) = {m_vals[3]}")
    print(f"M(5) = {m_vals[5]}")
    
    final_answer = f"{m_0},{m_vals[3]},{m_vals[5]}"
    print(f"\nFormatted answer: {final_answer}")
    
    # Final output as requested
    print(f"<<<{final_answer}>>>")


solve()
