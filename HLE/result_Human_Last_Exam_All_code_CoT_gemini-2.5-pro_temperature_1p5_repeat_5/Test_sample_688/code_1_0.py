import networkx as nx
from itertools import combinations
from math import factorial
from fractions import Fraction

def calculate_cn(n):
    """
    Calculates the prefactor c_n for a given number of nodes n.

    Args:
        n (int): The number of nodes in the diagram (must be >= 2).

    Returns:
        Fraction: The prefactor c_n as a fraction.
    """
    if n < 2:
        raise ValueError("n must be at least 2.")

    # A single edge is biconnected by convention in this context.
    if n == 2:
        # For n=2, the only biconnected graph is K_2 itself.
        # |E(K_2)| = 1. w_K2 = (-1)^(1-1) = 1.
        # c_2 = -(2-1)/2! * w_K2 = -1/2 * 1 = -1/2
        w_kn = 1
        print("For n = 2:")
        print("  The only biconnected graph is K_2 itself.")
        print(f"  w_2 = 1")
    else:
        # For n >= 3, biconnected requires at least 2 vertex-disjoint paths between any two vertices.
        k_n = nx.complete_graph(n)
        num_k_n_edges = k_n.number_of_edges()
        
        print(f"For n = {n}:")
        print(f"  Total edges in K_{n}: {num_k_n_edges}")
        
        w_kn = 0
        biconnected_counts = {}

        # Iterate over all possible subgraphs of K_n
        # A subgraph is defined by its edge set.
        all_edges = list(k_n.edges())
        for k in range(1, num_k_n_edges + 1):
            for edge_subset in combinations(all_edges, k):
                # Create the subgraph
                G = nx.Graph()
                G.add_nodes_from(range(n))
                G.add_edges_from(edge_subset)

                # Check for biconnectivity
                if nx.is_biconnected(G):
                    w_kn += (-1)**(num_k_n_edges - G.number_of_edges())
                    num_edges = G.number_of_edges()
                    biconnected_counts[num_edges] = biconnected_counts.get(num_edges, 0) + 1
        
        print("  Biconnected subgraphs found:")
        for edges, count in sorted(biconnected_counts.items()):
            print(f"    - {edges} edges: {count} graphs")
        
        print(f"  w_{n} = {w_kn}")


    # Calculate c_n = -(n-1)/n! * w_kn
    c_n_frac = -Fraction(n - 1, factorial(n)) * w_kn

    # Output the final equation
    # Using format to show numerator and denominator explicitly.
    num = (n - 1)
    den = factorial(n)
    w_val = w_kn
    
    print(f"\n  The prefactor c_{n} is calculated as:")
    print(f"  c_{n} = -({n}-1)/{n}! * w_{n}")
    print(f"  c_{n} = -({num})/{den} * ({w_val})")
    print(f"  c_{n} = {c_n_frac.numerator}/{c_n_frac.denominator}")
    print("-" * 20)
    return c_n_frac

if __name__ == '__main__':
    # Determine the prefactor for several values of n
    # For n > 6, the calculation becomes very slow.
    for n_val in range(2, 6):
        calculate_cn(n_val)
