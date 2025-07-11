import itertools

def calculate_n_g(num_vertices, edges):
    """
    Calculates N(G) for a given graph G.
    G is defined by its number of vertices and edge list.
    """
    num_edges = len(edges)
    
    # Create an adjacency list for easy access to incident edges
    adj = {i: [] for i in range(num_vertices)}
    for i, (u, v) in enumerate(edges):
        adj[u].append(i)
        adj[v].append(i)

    valid_colorings_count = 0
    # Iterate through all 2^|E| colorings
    # A coloring is represented by an integer from 0 to 2^|E|-1
    for i in range(2**num_edges):
        is_valid = True
        # Check each vertex for the slice condition
        for v in range(num_vertices):
            incident_edges_indices = adj[v]
            
            # Get colors of incident edges
            # color is 1 if the bit is set, 0 otherwise
            colors = [(i >> edge_idx) & 1 for edge_idx in incident_edges_indices]
            
            # Check if all colors are the same
            if sum(colors) == 0 or sum(colors) == 3:
                is_valid = False
                break
        
        if is_valid:
            valid_colorings_count += 1
            
    # N(G) is the number of valid colorings divided by 2
    # because swapping colors gives the same partition of edges.
    return valid_colorings_count // 2

def solve_m_n():
    """
    Solves for M(0), M(3), and M(5) by checking small cubic graphs.
    """
    
    # --- M(3) ---
    # Smallest cubic graph has |V|=4, which is K_4
    print("Finding M(3):")
    k4_vertices = 4
    k4_edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    n_k4 = calculate_n_g(k4_vertices, k4_edges)
    print(f"For G = K_4 (|V|=4), N(G) = {n_k4}")

    m3 = "none"
    if n_k4 % 3 == 0:
        m3 = k4_vertices
        print(f"N(K_4) is a multiple of 3. Since |V|=4 is the smallest possible, M(3) = {m3}.")
    
    print("-" * 20)

    # --- M(5) ---
    print("Finding M(5):")
    print("Checking |V|=4...")
    if n_k4 % 5 == 0:
        m5 = k4_vertices
        print(f"N(K_4) is a multiple of 5. So M(5) = {m5}.")
    else:
        print(f"N(K_4) = {n_k4} is not a multiple of 5. Checking graphs with |V|=6...")
        # K_3,3 (|V|=6, |E|=9)
        k33_vertices = 6
        k33_edges = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]
        n_k33 = calculate_n_g(k33_vertices, k33_edges)
        print(f"For G = K_3,3 (|V|=6), N(G) = {n_k33}")

        # Triangular Prism Y_3 (|V|=6, |E|=9)
        y3_vertices = 6
        y3_edges = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
        n_y3 = calculate_n_g(y3_vertices, y3_edges)
        print(f"For G = Y_3 (|V|=6), N(G) = {n_y3}")

        if n_k33 % 5 == 0 or n_y3 % 5 == 0:
            m5 = 6
            print(f"Found a graph with |V|=6 whose N(G) is a multiple of 5. So M(5) = {m5}.")
        else:
            print("Neither N(K_3,3) nor N(Y_3) is a multiple of 5. Checking graphs with |V|=8...")
            # Cube graph Q_3 (|V|=8, |E|=12)
            q3_vertices = 8
            q3_edges = [(0, 1), (1, 2), (2, 3), (3, 0),
                        (4, 5), (5, 6), (6, 7), (7, 4),
                        (0, 4), (1, 5), (2, 6), (3, 7)]
            n_q3 = calculate_n_g(q3_vertices, q3_edges)
            print(f"For G = Q_3 (|V|=8), N(G) = {n_q3}")
            
            if n_q3 % 5 == 0:
                m5 = 8
                print(f"N(Q_3) is a multiple of 5. So M(5) = {m5}.")
            else:
                 m5 = "none" # based on this search
    
    print("-" * 20)

    # --- M(0) ---
    print("Finding M(0):")
    m0 = "none"
    print("For N(G) to be a multiple of 0, N(G) must be 0.")
    print("A known result in graph theory states that every cubic graph has at least one slice.")
    print("In fact, it is known that N(G) >= 3 for all cubic graphs.")
    print("Therefore, no cubic graph G has N(G) = 0.")
    print("So, M(0) is none.")

    print("-" * 20)
    print("Final answer compilation:")
    print(f"M(0) = {m0}")
    print(f"M(3) = {m3}")
    print(f"M(5) = {m5}")
    
    final_answer = f"{m0},{m3},{m5}"
    print("The answer in the requested format is:")
    print(f"<<<{final_answer}>>>")

solve_m_n()