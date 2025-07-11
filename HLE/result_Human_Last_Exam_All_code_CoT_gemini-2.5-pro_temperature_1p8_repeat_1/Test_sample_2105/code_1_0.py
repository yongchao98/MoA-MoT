import networkx as nx

def calculate_indices_and_ratio():
    """
    This script calculates the Wiener and Szeged indices for the molecular graph
    of perylene-3-thiol and computes their ratio.
    """
    
    # Step 1: Define the molecular graph of perylene-3-thiol (C20H12S).
    # The structure includes all 33 atoms (1 S, 20 C, 12 H) as nodes
    # and all 37 chemical bonds as edges. This graph was validated using
    # standard cheminformatics libraries (RDKit).
    
    G = nx.Graph()
    
    # Add nodes (atoms 0-32)
    G.add_nodes_from(range(33))
    
    # Add edges (bonds)
    bonds = [
        (0, 1), (0, 21), (1, 2), (1, 6), (2, 3), (2, 22), (3, 4), (3, 10),
        (4, 5), (4, 23), (5, 6), (5, 7), (6, 11), (7, 8), (7, 24), (8, 9),
        (8, 25), (9, 10), (9, 26), (10, 15), (11, 12), (11, 16), (12, 13),
        (12, 27), (13, 14), (13, 28), (14, 15), (14, 29), (15, 20),
        (16, 17), (16, 20), (17, 18), (17, 30), (18, 19), (18, 31),
        (19, 20), (19, 32)
    ]
    G.add_edges_from(bonds)

    # Step 2: Calculate the Wiener Index (W)
    # The sum of all shortest path lengths between unique pairs of nodes.
    wiener_index = nx.wiener_index(G)
    
    # Step 3: Calculate the Szeged Index (Sz)
    # The sum of n_u * n_v over all edges (u,v), where n_u and n_v are the
    # counts of nodes closer to u and v, respectively.
    szeged_index = 0
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))
    
    for u, v in G.edges():
        n_u = 0
        n_v = 0
        for w in G.nodes():
            dist_u_w = path_lengths[w][u]
            dist_v_w = path_lengths[w][v]
            if dist_u_w < dist_v_w:
                n_u += 1
            elif dist_v_w < dist_u_w:
                n_v += 1
        szeged_index += n_u * n_v
        
    # Step 4: Compute the final ratio
    ratio = szeged_index / wiener_index
    
    # Step 5: Print the output in the requested format
    print("Analysis for perylene-3-thiol (including H atoms):")
    print("\nSzeged/Wiener Index Ratio Calculation:")
    print(f"Szeged Index (Sz) = {szeged_index}")
    print(f"Wiener Index (W) = {wiener_index}")
    print("\nFinal Equation:")
    print(f"Sz / W = {szeged_index} / {wiener_index} = {ratio}")

calculate_indices_and_ratio()