import networkx as nx

def create_lattice_graph():
    """
    Creates the Lattice graph L_2(4), also known as the 4x4 Rook's graph.
    This is one of the two SRGs with parameters (16, 6, 2, 2).
    Vertices are cells on a 4x4 grid, adjacent if they share a row or column.
    """
    G = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G.add_nodes_from(nodes)

    # Add edges for vertices in the same row
    for r in range(4):
        for c1 in range(4):
            for c2 in range(c1 + 1, 4):
                G.add_edge((r, c1), (r, c2))

    # Add edges for vertices in the same column
    for c in range(4):
        for r1 in range(4):
            for r2 in range(r1 + 1, 4):
                G.add_edge((r1, c), (r2, c))

    return G

def create_shrikhande_graph():
    """
    Creates the Shrikhande graph.
    This is the other SRG with parameters (16, 6, 2, 2).
    Vertices are on a 4x4 grid with wrap-around (torus) connections.
    A vertex (i,j) is adjacent to its six neighbors:
    (i +/- 1, j), (i, j +/- 1), (i+1, j+1), and (i-1, j-1), all modulo 4.
    """
    G = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G.add_nodes_from(nodes)
    
    # Define the six offsets for neighbors on the Z_4 x Z_4 grid
    neighbor_offsets = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]
    
    for i in range(4):
        for j in range(4):
            v1 = (i, j)
            for dx, dy in neighbor_offsets:
                # Calculate neighbor coordinates with wrap-around
                v2 = ((i + dx) % 4, (j + dy) % 4)
                G.add_edge(v1, v2)
    return G

def count_simple_cycles_of_length_k(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    This is done via a recursive backtracking search (DFS) starting from each node.
    """
    # Using the built-in simple_cycles generator from networkx is a robust way to do this.
    # It finds all simple cycles, and we filter for those of the desired length.
    cycle_generator = nx.simple_cycles(G)
    count = 0
    for cycle in cycle_generator:
        if len(cycle) == k:
            count += 1
    return count

# --- Main Execution ---
print("This script checks if two SRGs with the same parameters can have a different number of 5-cycles.")
print("The chosen parameters are (n, d, lambda, mu) = (16, 6, 2, 2).\n")

# Create the two graphs
G_lattice = create_lattice_graph()
G_shrikhande = create_shrikhande_graph()

print("Graph 1: The Lattice graph L_2(4)")
print("Graph 2: The Shrikhande graph\n")

# Count the 5-cycles in each graph
# This operation can be computationally intensive and may take a moment.
print("Counting 5-cycles in the Lattice graph...")
c5_lattice = count_simple_cycles_of_length_k(G_lattice, 5)

print("Counting 5-cycles in the Shrikhande graph...")
c5_shrikhande = count_simple_cycles_of_length_k(G_shrikhande, 5)

print("\n--- Results ---")
print(f"Number of 5-cycles in the Lattice graph: {c5_lattice}")
print(f"Number of 5-cycles in the Shrikhande graph: {c5_shrikhande}")

if c5_lattice != c5_shrikhande:
    print("\nThe number of 5-cycles is different, so such a combination of parameters and graphs exists.")
else:
    print("\nThe number of 5-cycles is the same for this pair.")

<<<Yes>>>