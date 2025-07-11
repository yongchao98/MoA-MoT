import networkx as nx

def calculate_indices():
    """
    Calculates and prints the Wiener and Szeged indices, and their ratio,
    for the major reduction product of di(perylene-3-yl) disulfide, which is perylene-3-thiol.
    """
    # Step 1: Create the carbon skeleton of perylene using the built-in generator.
    # This creates a graph with 20 nodes representing the carbon atoms.
    G = nx.perylenegraph()

    # Step 2: Add hydrogens and the thiol group (-SH) to the skeleton.
    # In the perylene skeleton, carbons that bind to hydrogens have a degree of 2.
    # Perylene (C20H12) has 12 such carbons.
    degree_2_carbons = [node for node, degree in G.degree() if degree == 2]

    # Define indices for the new atoms to be added.
    # Carbons are 0-19. Sulfur will be 20. Hydrogens will be 21-32.
    s_atom_idx = 20
    h_thiol_idx = 21

    # In perylene-3-thiol, one C-H bond is replaced by a C-SH group.
    # Due to the high symmetry of perylene, we can substitute any of the
    # 12 equivalent positions. We pick the first one.
    c_thiol_idx = degree_2_carbons[0]
    G.add_edge(c_thiol_idx, s_atom_idx)
    G.add_edge(s_atom_idx, h_thiol_idx)

    # Add the remaining 11 hydrogen atoms to the other degree-2 carbons.
    h_carbons = degree_2_carbons[1:]
    h_atom_indices = range(22, 33)  # Indices for the 11 H atoms on the ring
    for c_idx, h_idx in zip(h_carbons, h_atom_indices):
        G.add_edge(c_idx, h_idx)

    # The graph now represents perylene-3-thiol with 33 atoms.

    # Step 3: Calculate the Wiener Index (W).
    # W = 1/2 * sum(d(u, v)) for all pairs of nodes u, v.
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))
    wiener_index = sum(sum(lengths.values()) for lengths in path_lengths.values()) / 2

    # Step 4: Calculate the Szeged Index (Sz).
    # Sz = sum over all edges e=(u,v) of n_u(e) * n_v(e).
    szeged_index = 0
    for u, v in G.edges():
        n_u = 0
        n_v = 0
        path_lengths_u = path_lengths[u]
        path_lengths_v = path_lengths[v]
        for w in G.nodes():
            if path_lengths_u[w] < path_lengths_v[w]:
                n_u += 1
            elif path_lengths_v[w] < path_lengths_u[w]:
                n_v += 1
        szeged_index += n_u * n_v

    # Step 5: Print the results and the final equation.
    # The indices must be integers.
    W = int(wiener_index)
    Sz = int(szeged_index)
    ratio = Sz / W if W != 0 else 0
    
    print(f"The major reduction product is perylene-3-thiol (C20H11SH).")
    print(f"The molecular graph has {G.number_of_nodes()} atoms (nodes) and {G.number_of_edges()} bonds (edges).")
    print("-" * 30)
    print(f"Wiener Index (W): {W}")
    print(f"Szeged Index (Sz): {Sz}")
    print(f"Szeged/Wiener Ratio calculation:")
    print(f"{Sz} / {W} = {ratio}")


if __name__ == '__main__':
    calculate_indices()

<<<3.7041564792176245>>>