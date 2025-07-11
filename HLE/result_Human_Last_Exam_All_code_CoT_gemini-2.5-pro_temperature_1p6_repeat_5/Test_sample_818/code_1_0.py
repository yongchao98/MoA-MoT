import itertools

def solve():
    """
    This function calculates N(K_4) and then determines M(0), M(3), and M(5).
    """

    # K_4 has 4 vertices and 6 edges.
    num_vertices = 4
    vertices = list(range(num_vertices))
    edges = list(itertools.combinations(vertices, 2))
    num_edges = len(edges)

    # Adjacency list representation for K_4
    # For each vertex, store the indices of its incident edges
    incident_edges = [[] for _ in range(num_vertices)]
    for i, edge in enumerate(edges):
        u, v = edge
        incident_edges[u].append(i)
        incident_edges[v].append(i)

    slice_colorings_count = 0
    sc_slice_colorings_count = 0 # Self-complementary slice colorings

    # Iterate through all 2^6 = 64 possible edge colorings
    # A coloring is represented by an integer from 0 to 63
    # Bit `j` of the integer represents the color of edge `j` (0 or 1)
    for i in range(2**num_edges):
        is_slice = True
        # Check each vertex for the slice condition
        for v_idx in range(num_vertices):
            
            # Get colors of incident edges
            colors = []
            for edge_idx in incident_edges[v_idx]:
                # Check j-th bit of i
                if (i >> edge_idx) & 1:
                    colors.append(1)
                else:
                    colors.append(0)
            
            # If all incident edges have the same color, it's not a slice
            if all(c == colors[0] for c in colors):
                is_slice = False
                break
        
        if is_slice:
            slice_colorings_count += 1
            
            # Check if the coloring is self-complementary
            # Number of edges with color 1
            num_color_1 = bin(i).count('1')
            if num_color_1 == num_edges / 2:
                sc_slice_colorings_count += 1

    # The number of unique slices N(G) is N_sc + (S - N_sc) / 2
    # S = slice_colorings_count
    # N_sc = sc_slice_colorings_count
    n_k4 = sc_slice_colorings_count + (slice_colorings_count - sc_slice_colorings_count) // 2

    # Now, we determine M(0), M(3), M(5)
    
    # M(0): As argued in the plan, every cubic graph has a slice, so N(G) >= 1.
    # N(G) cannot be 0, so M(0) is none.
    m_0 = "none"

    # M(3) and M(5): The smallest cubic graph is K_4 (m=4).
    # We calculated N(K_4).
    m_3 = "none"
    if n_k4 % 3 == 0:
        m_3 = 4

    m_5 = "none"
    if n_k4 % 5 == 0:
        m_5 = 4
        
    # We calculated the following values inside the code:
    # slice_colorings_count = 18
    # sc_slice_colorings_count = 12
    # n_k4 = 12 + (18 - 12) // 2 = 15

    # This part of the code prints the reasoning and the final answer.
    print(f"For K_4 (m=4):")
    print(f"  - Total slice colorings (S): {slice_colorings_count}")
    print(f"  - Self-complementary slice colorings (N_sc): {sc_slice_colorings_count}")
    print(f"  - Number of slices N(K_4) = N_sc + (S - N_sc)/2 = {n_k4}")
    print("\nDetermining M(n):")
    print(f"M(0): N(G) can never be 0 for any cubic graph G. Thus M(0) is 'none'.")
    print(f"M(3): For m=4, N(K_4)={n_k4}. Since {n_k4} is a multiple of 3, and m=4 is the smallest possible, M(3) = 4.")
    print(f"M(5): For m=4, N(K_4)={n_k4}. Since {n_k4} is a multiple of 5, and m=4 is the smallest possible, M(5) = 4.")

    final_answer = f"{m_0},{m_3},{m_5}"
    print(f"\nFinal answer in the required format: {final_answer}")
    print(f'<<<{final_answer}>>>')

solve()