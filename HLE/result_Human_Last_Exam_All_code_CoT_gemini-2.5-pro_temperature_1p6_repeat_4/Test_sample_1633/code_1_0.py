import networkx as nx
import random
import numpy as np

def demonstrate_clustering_impact():
    """
    Demonstrates the impact on the clustering coefficient of removing
    a lattice edge vs. a shortcut edge in a Watts-Strogatz graph.
    """
    n = 1000  # Number of nodes
    k = 6     # Initial degree (k_0)
    p = 0.2   # Rewiring probability (beta)

    # Generate the Watts-Strogatz graph
    # networkx generates the lattice and then rewires, so we can find
    # which edges were rewired vs. which are original lattice edges.
    G = nx.watts_strogatz_graph(n, k, p, seed=42)

    # Initial clustering coefficient
    c_initial = nx.average_clustering(G)
    print(f"Initial Watts-Strogatz Graph (n={n}, k={k}, p={p})")
    print(f"Initial Average Clustering Coefficient: {c_initial:.4f}\n")

    # The WS model in networkx starts with a k-regular ring lattice.
    # We can identify lattice edges as those connecting nearby nodes.
    # An edge (u, v) is a lattice edge if |u-v| mod n is small.
    # Shortcuts are edges where |u-v| mod n is large.
    
    shortcuts = []
    lattice_edges = []
    # k/2 is the number of neighbors on each side
    # Any distance > k/2 is not an immediate lattice neighbor
    for u, v in G.edges():
        dist = min(abs(u - v), n - abs(u - v))
        if dist > k / 2:
            shortcuts.append((u, v))
        else:
            lattice_edges.append((u, v))

    print(f"Found {len(lattice_edges)} lattice edges and {len(shortcuts)} shortcut edges.")
    
    if not lattice_edges or not shortcuts:
        print("Could not find both edge types to compare. Exiting.")
        return

    # --- Case 1: Remove a random lattice edge ---
    edge_to_remove_lattice = random.choice(lattice_edges)
    G_removed_lattice = G.copy()
    G_removed_lattice.remove_edge(*edge_to_remove_lattice)
    c_after_lattice_removal = nx.average_clustering(G_removed_lattice)
    
    print(f"\n--- Removing a 'lattice' edge {edge_to_remove_lattice} ---")
    print(f"Clustering Coefficient after removal: {c_after_lattice_removal:.4f}")
    print(f"Drop in C: {c_initial - c_after_lattice_removal:.6f}")

    # --- Case 2: Remove a random shortcut edge ---
    edge_to_remove_shortcut = random.choice(shortcuts)
    G_removed_shortcut = G.copy()
    G_removed_shortcut.remove_edge(*edge_to_remove_shortcut)
    c_after_shortcut_removal = nx.average_clustering(G_removed_shortcut)

    print(f"\n--- Removing a 'shortcut' edge {edge_to_remove_shortcut} ---")
    print(f"Clustering Coefficient after removal: {c_after_shortcut_removal:.4f}")
    print(f"Drop in C: {c_initial - c_after_shortcut_removal:.6f}")

    # The drop for removing a lattice edge is typically much higher
    # because that edge is more likely to be part of multiple triangles.
    # This supports the argument that Option H is unlikely to be correct
    # as it would violate the clustering constraint.

if __name__ == '__main__':
    demonstrate_clustering_impact()
