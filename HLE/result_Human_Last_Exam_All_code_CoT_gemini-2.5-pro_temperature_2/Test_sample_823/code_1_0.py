import networkx as nx
import numpy as np

def analyze_grid_graphs(k_values):
    """
    Analyzes k-by-k grid graphs to demonstrate a counterexample for option C.
    A family of expander graphs must have its algebraic connectivity bounded away from zero.
    This function will show that for the family of grid graphs, the algebraic connectivity
    approaches zero as k increases, while treewidth is unbounded and degree is bounded.
    """
    print("Analyzing properties of k-by-k grid graphs:")
    print("-" * 70)
    print(f"{'k':<5} | {'Num Vertices':<15} | {'Max Degree':<12} | {'Treewidth':<10} | {'Algebraic Connectivity':<25}")
    print("-" * 70)

    for k in k_values:
        # Create a k-by-k grid graph
        G = nx.grid_2d_graph(k, k)

        # Number of vertices
        num_vertices = G.number_of_nodes()

        # Maximum degree in a grid graph is at most 4
        # We find the max degree for confirmation
        degrees = [d for n, d in G.degree()]
        max_degree = max(degrees) if degrees else 0

        # The treewidth of a k-by-k grid is k
        treewidth = k

        # Calculate the algebraic connectivity (Fiedler value)
        # For non-connected graphs, it is 0. Grids are connected.
        # This requires scipy to be installed.
        try:
            alg_conn = nx.algebraic_connectivity(G)
        except ImportError:
            print("Scipy is required to calculate algebraic connectivity.")
            return
        except np.linalg.LinAlgError:
            # This may happen for very small k, e.g., k=1
            alg_conn = 0

        print(f"{k:<5} | {num_vertices:<15} | {max_degree:<12} | {treewidth:<10} | {alg_conn:<25.5f}")
    print("-" * 70)
    print("\nConclusion: The family of grid graphs has bounded degree (4) and unbounded")
    print("treewidth (k), but its algebraic connectivity approaches 0. Therefore, it is")
    print("not an expander family, which falsifies option C.")


if __name__ == '__main__':
    # Define a series of k values to test
    k_series = [5, 10, 15, 20, 25]
    analyze_grid_graphs(k_series)