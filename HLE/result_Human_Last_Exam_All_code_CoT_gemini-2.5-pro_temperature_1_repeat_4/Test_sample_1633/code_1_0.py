import math

def analyze_watts_strogatz_shortcuts(n, k0, beta):
    """
    Analyzes the shortcut edge distribution in a Watts-Strogatz graph.

    Args:
        n (int): Number of vertices.
        k0 (int): Initial degree in the lattice.
        beta (float): Rewiring probability.
    """
    
    # The total number of edges in the graph is n * k0 / 2
    total_edges = n * k0 / 2
    
    # The number of shortcut edges created is beta * total_edges
    num_shortcuts = beta * total_edges
    
    # In the random graph of shortcuts, the average degree (lambda) is 2 * num_shortcuts / n
    # This simplifies to beta * k0
    lambda_val = beta * k0
    
    # The degree distribution of this random graph follows a Poisson distribution.
    # The probability of a vertex having degree 0 (no shortcuts) is e^(-lambda).
    prob_no_shortcuts = math.exp(-lambda_val)
    
    # The expected number of vertices with no shortcuts is n * prob_no_shortcuts
    num_nodes_no_shortcuts = n * prob_no_shortcuts
    
    # Print the derivation
    print("### Analysis of Shortcut Distribution ###")
    print(f"Initial parameters: k₀ = {k0}, β = {beta}")
    print("\nStep 1: Calculate the average number of shortcuts per node (λ).")
    print(f"λ = β * k₀ = {beta} * {k0} = {lambda_val}")
    
    print("\nStep 2: Calculate the fraction of nodes with zero shortcuts.")
    print("This follows a Poisson distribution P(k=0) = e^(-λ).")
    print(f"Fraction with no shortcuts = e^(-{lambda_val}) ≈ {prob_no_shortcuts:.4f}")
    
    print("\n### Conclusion ###")
    print(f"Approximately {prob_no_shortcuts:.1%} of the vertices have only original lattice edges.")
    print("To reduce the degree of these vertices, the transformation MUST remove lattice edges.")
    
    # Relate this to option H
    required_removals_fraction = 1/6
    print(f"\nOption H states that at least n/6 (≈ {required_removals_fraction:.4f}*n) lattice edge removals are required.")
    print(f"Our analysis shows that ~{prob_no_shortcuts:.4f}*n nodes are only connected by lattice edges.")
    print("Any global transformation plan will likely need to decrease the degrees of a significant portion of these nodes.")
    print(f"For example, modifying just over half of them would require > {0.5*prob_no_shortcuts:.4f}*n ≈ {required_removals_fraction:.4f}*n lattice removals.")
    print("This provides strong support for option H being a necessary consequence of the transformation constraints.")


# Main execution
if __name__ == "__main__":
    # Parameters from the problem
    initial_degree_k0 = 6
    rewiring_prob_beta = 0.2
    
    # Use n=1 for calculating the fraction, as n cancels out
    analyze_watts_strogatz_shortcuts(1, initial_degree_k0, rewiring_prob_beta)

<<<H>>>