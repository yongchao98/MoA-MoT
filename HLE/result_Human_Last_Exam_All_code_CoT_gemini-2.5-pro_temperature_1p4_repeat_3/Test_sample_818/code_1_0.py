import itertools

def count_slices(name, edges, num_vertices):
    """
    Computes the number of slices N(G) for a given graph G.
    A slice is a partition of edges into two classes where each vertex is incident to at least one edge in each class.
    This is equivalent to half the number of 2-edge-colorings where no vertex is monochromatic.
    """
    num_edges = len(edges)
    adj = {i: [] for i in range(num_vertices)}
    for i, (u, v) in enumerate(edges):
        adj[u].append(i)
        adj[v].append(i)

    valid_colorings = 0
    # Iterate through all 2^|E| possible edge colorings.
    # A coloring is represented by an integer, where the k-th bit is the color of edge k.
    for i in range(2**num_edges):
        is_valid_coloring = True
        # Check the monochromatic condition at each vertex
        for v in range(num_vertices):
            # All graphs considered are cubic, so degree is 3.
            edge_indices = adj[v]
            c1 = (i >> edge_indices[0]) & 1
            c2 = (i >> edge_indices[1]) & 1
            c3 = (i >> edge_indices[2]) & 1
            if c1 == c2 and c2 == c3:
                is_valid_coloring = False
                break # This coloring is invalid, move to the next.
        
        if is_valid_coloring:
            valid_colorings += 1
            
    # The number of slices is half the number of valid colorings, since swapping the two colors
    # results in a valid coloring that corresponds to the same partition.
    num_slices = valid_colorings // 2
    print(f"Graph: {name} ({num_vertices} vertices). N({name}) = {num_slices}")
    return num_slices

def solve():
    """
    Finds M(0), M(3), and M(5).
    """
    # M(0)
    # A multiple of 0 must be 0. We need N(G) = 0.
    # However, a theorem by A. Holroyd states that every finite cubic graph has a slice, so N(G) >= 1.
    # Therefore, no such graph exists.
    m0 = "none"
    print("M(0): Based on the theorem that all cubic graphs have at least one slice, N(G) is never 0. So M(0) is none.")
    
    # M(3)
    # Start searching from the smallest number of vertices for a cubic graph, m=4.
    print("\nSearching for M(3)...")
    
    # m = 4: The complete graph K4 is the only cubic graph.
    k4_edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    n_k4 = count_slices("K_4", k4_edges, 4)
    if n_k4 % 3 == 0:
        m3 = 4
        print(f"N(K_4)={n_k4} is a multiple of 3. Smallest m is 4. M(3) = 4.")
    else:
        # This part of the code won't be reached, but is included for completeness.
        m3 = "not found yet"

    # M(5)
    print("\nSearching for M(5)...")
    # m = 4: K4
    if n_k4 % 5 == 0:
        m5 = 4
        print(f"N(K_4)={n_k4} is a multiple of 5. Smallest m is 4. M(5) = 4.")
    else:
        print(f"N(K_4)={n_k4} is not a multiple of 5. Checking m=6.")
        # m = 6: Prism graph and K_3,3
        prism_edges = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
        n_prism = count_slices("Prism(C3xK2)", prism_edges, 6)

        k33_edges = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]
        n_k33 = count_slices("K_3,3", k33_edges, 6)
        
        if n_prism % 5 == 0 or n_k33 % 5 == 0:
            m5 = 6
            print(f"Found a graph on 6 vertices whose N(G) is a multiple of 5. M(5) = 6.")
        else:
            print(f"N(Prism)={n_prism} and N(K_3,3)={n_k33} are not multiples of 5. Checking m=8.")
            # m = 8: There are 5 cubic graphs. We only need to check until one works. The Cube graph is one of them.
            cube_edges = [(0,1), (1,2), (2,3), (3,0), (4,5), (5,6), (6,7), (7,4), (0,4), (1,5), (2,6), (3,7)]
            n_cube = count_slices("Cube(Q3)", cube_edges, 8)
            if n_cube % 5 == 0:
                m5 = 8
                print(f"N(Cube)={n_cube} is a multiple of 5. Smallest m is 8. M(5) = 8.")
            else:
                 m5 = "not found yet"
    
    print(f"\nFinal Answer: M(0), M(3), M(5) = {m0},{m3},{m5}")
    return f"{m0},{m3},{m5}"

if __name__ == '__main__':
    final_answer = solve()
    # The final answer format required by the prompt
    # print(f"<<<{final_answer}>>>")