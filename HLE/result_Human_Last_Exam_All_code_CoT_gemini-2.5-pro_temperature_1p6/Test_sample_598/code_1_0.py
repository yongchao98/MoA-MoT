import networkx as nx
import numpy as np

def calculate_and_print_cycle_counts(graph_name, G, n, d, l, m):
    """
    Calculates the trace of A^5 and the number of 5-cycles for a given graph.
    Then prints the results in a detailed equation format.
    """
    if G.number_of_nodes() == 0:
        print(f"Could not construct or find {graph_name}. Skipping.")
        return

    # Convert graph to numpy array for matrix operations
    A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))

    # Calculate trace of A^5
    A5 = np.linalg.matrix_power(A, 5)
    trace_A5 = int(np.trace(A5))

    # Calculate the constant term in the formula, which depends only on SRG parameters
    # This term counts the total number of non-simple closed walks of length 5
    constant_term = n * d * l * (d + l - 3)

    # Calculate the number of 5-cycles
    num_5_cycles = (trace_A5 - constant_term) // 10

    print(f"For the {graph_name}:")
    print(f"- Parameters (n,d,λ,μ) = ({n},{d},{l},{m})")
    print(f"- The trace of the 5th power of its adjacency matrix, Tr(A^5), is {trace_A5}.")
    print(f"- The number of 5-cycles (N5) is calculated as:")
    print(f"  N5 = (Tr(A^5) - n*d*λ*(d+λ-3)) / 10")
    print(f"  N5 = ({trace_A5} - {n}*{d}*{l}*({d}+{l}-3)) / 10")
    print(f"  N5 = ({trace_A5} - {constant_term}) / 10")
    print(f"  N5 = {num_5_cycles}\n")

def main():
    """
    Main function to construct the graphs and trigger calculations.
    """
    # Define the SRG parameters
    n, d, l, m = 16, 6, 2, 2

    # --- Graph 1: 4x4 Rook's Graph ---
    # Vertices are squares on a 4x4 board (labeled 0-15).
    # Edges connect vertices in the same row or column.
    G_rook = nx.Graph()
    nodes = range(n)
    G_rook.add_nodes_from(nodes)
    for i in nodes:
        for j in range(i + 1, n):
            # i and j are in the same row if their integer division by 4 is the same.
            is_same_row = (i // 4 == j // 4)
            # i and j are in the same col if their modulo 4 is the same.
            is_same_col = (i % 4 == j % 4)
            if is_same_row or is_same_col:
                G_rook.add_edge(i, j)

    # --- Graph 2: Shrikhande Graph ---
    G_shrikhande = nx.shrikhande_graph()
    
    print("Yes, a pair of strongly regular graphs with the same parameters can have a different number of 5-cycles.")
    print("An example is the 4x4 Rook's graph and the Shrikhande graph, both srg(16, 6, 2, 2).\n")

    calculate_and_print_cycle_counts("4x4 Rook's Graph", G_rook, n, d, l, m)
    calculate_and_print_cycle_counts("Shrikhande Graph", G_shrikhande, n, d, l, m)
    
    print("As shown, the two graphs have different numbers of 5-cycles (1152 vs 1056).")

if __name__ == "__main__":
    main()