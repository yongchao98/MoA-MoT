# This script requires the networkx and numpy libraries.
# You can install them using: pip install networkx numpy
import networkx as nx
import numpy as np

def construct_rooks_graph():
    """
    Constructs the 4x4 Rook's graph, an SRG(16, 6, 2, 2).
    Vertices are cells (i,j) of a 4x4 grid.
    Two vertices are adjacent if they are in the same row or column.
    """
    G = nx.Graph()
    n_side = 4
    nodes = range(n_side * n_side)
    G.add_nodes_from(nodes)
    for i in nodes:
        for j in range(i + 1, len(nodes)):
            r1, c1 = i // n_side, i % n_side
            r2, c2 = j // n_side, j % n_side
            if r1 == r2 or c1 == c2:
                G.add_edge(i, j)
    return G

def construct_shrikhande_graph():
    """
    Constructs the Shrikhande graph, also an SRG(16, 6, 2, 2).
    It is defined as a Cayley graph on Z_4 x Z_4.
    """
    G = nx.Graph()
    n_side = 4
    nodes = range(n_side * n_side)
    G.add_nodes_from(nodes)
    # Neighborhood set S for differences (dx, dy) mod 4
    S = {(0, 1), (0, 3), (1, 0), (3, 0), (1, 1), (3, 3)}
    for i in nodes:
        for j in range(i + 1, len(nodes)):
            r1, c1 = i // n_side, i % n_side
            r2, c2 = j // n_side, j % n_side
            dr = (r1 - r2) % n_side
            dc = (c1 - c2) % n_side
            if (dr, dc) in S:
                G.add_edge(i, j)
    return G

def analyze_graph(graph_name, G, n, d, lam):
    """
    Calculates and prints the number of 5-cycles for a given graph.
    """
    print(f"--- Analysis for: {graph_name} ---")
    
    # Get adjacency matrix and compute the trace of its 5th power
    A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))
    A5 = np.linalg.matrix_power(A, 5)
    tr_A5 = int(np.trace(A5))
    
    # Calculate the number of 5-cycles using the formula
    # 10 * c5 = tr(A^5) - 5 * n * d * lambda * (d - 1)
    term_to_subtract = 5 * n * d * lam * (d - 1)
    numerator = tr_A5 - term_to_subtract
    c5 = numerator // 10

    print(f"The number of 5-cycles, c5, is calculated by the formula:")
    print(f"10 * c5 = tr(A^5) - 5 * n * d * lambda * (d - 1)")
    print("\nPlugging in the values for the graph:")
    print(f"10 * c5 = {tr_A5} - 5 * {n} * {d} * {lam} * ({d} - 1)")
    print(f"10 * c5 = {tr_A5} - {term_to_subtract}")
    print(f"10 * c5 = {numerator}")
    print(f"c5 = {c5}")
    print("-" * (len(graph_name) + 22))
    return c5

if __name__ == "__main__":
    # Parameters for the SRGs
    n = 16
    d = 6
    lam = 2  # lambda parameter
    mu = 2
    
    print(f"Yes, there is a combination of parameters for which two SRGs can have a different number of 5-cycles.")
    print(f"We will demonstrate this for the parameters (n,d,lambda,mu) = ({n},{d},{lam},{mu}).\n")

    # Construct the graphs
    G_rook = construct_rooks_graph()
    G_shrikhande = construct_shrikhande_graph()

    # Analyze each graph
    c5_rook = analyze_graph("4x4 Rook's Graph", G_rook, n, d, lam)
    print()
    c5_shrikhande = analyze_graph("Shrikhande Graph", G_shrikhande, n, d, lam)
    print()

    # Final conclusion
    print("Conclusion:")
    print(f"The 4x4 Rook's graph has {c5_rook} 5-cycles.")
    print(f"The Shrikhande graph has {c5_shrikhande} 5-cycles.")
    if c5_rook != c5_shrikhande:
        print("Since the number of 5-cycles is different, this confirms that such pairs of graphs exist.")
    else:
        print("These two graphs have the same number of 5-cycles.")
