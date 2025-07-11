def potts_interaction_count(config, graph_edges):
    """
    Calculates the interaction term S(config) for a Potts model configuration.
    S(config) = sum_{<u,v> in E} I(config[u] == config[v])
    
    Args:
        config (dict): A dictionary mapping vertex indices to spin states.
        graph_edges (list of tuples): A list where each tuple represents an edge.
    
    Returns:
        int: The sum of indicator functions, representing the number of monochromatic edges.
    """
    count = 0
    for u, v in graph_edges:
        if config[u] == config[v]:
            count += 1
    return count

def meet(xi, eta):
    """Computes the meet (coordinate-wise min) of two configurations."""
    keys = xi.keys()
    result = {v: min(xi[v], eta[v]) for v in keys}
    return result

def join(xi, eta):
    """Computes the join (coordinate-wise max) of two configurations."""
    keys = xi.keys()
    result = {v: max(xi[v], eta[v]) for v in keys}
    return result

def main():
    """
    This script demonstrates that the positive correlations property (FKG inequality)
    does not hold for the Potts model for any graph with at least one edge if q >= 3.
    It does this by showing a counterexample to the equivalent FKG lattice condition.
    The failure for a graph with maximum degree d=1 implies the answer to the question must be d=0.
    """
    print("Testing the FKG lattice condition for the Potts model.")
    print("The condition is: S(xi_meet_eta) + S(xi_join_eta) >= S(xi) + S(eta)")
    print("where S(config) is the number of edges with same-spin vertices.\n")
    
    # We choose the simplest graph with an edge: K_2. This graph has max degree d=1.
    graph_edges = [(0, 1)]
    print(f"Graph: A single edge between vertices 0 and 1. Maximum degree d=1.")
    
    # We choose q=3, so states can be {1, 2, 3}.
    print("Number of states q=3. The states are {1, 2, 3}.")

    # We choose specific configurations xi and eta that form a counterexample.
    xi = {0: 2, 1: 2}
    eta = {0: 1, 1: 3}
    print(f"Configuration xi: {xi}")
    print(f"Configuration eta: {eta}\n")

    # Compute the meet and join configurations
    xi_meet_eta = meet(xi, eta)
    xi_join_eta = join(xi, eta)
    print(f"Meet configuration (xi ^ eta): {xi_meet_eta}")
    print(f"Join configuration (xi v eta): {xi_join_eta}\n")
    
    # Calculate the interaction counts (S values)
    s_xi = potts_interaction_count(xi, graph_edges)
    s_eta = potts_interaction_count(eta, graph_edges)
    s_meet = potts_interaction_count(xi_meet_eta, graph_edges)
    s_join = potts_interaction_count(xi_join_eta, graph_edges)
    
    # Evaluate both sides of the inequality
    lhs_val = s_meet + s_join
    rhs_val = s_xi + s_eta
    
    print("Evaluating the two sides of the inequality 'S(xi^eta) + S(xiveta) >= S(xi) + S(eta)':")
    print(f"LHS = S({xi_meet_eta}) + S({xi_join_eta}) = {s_meet} + {s_join} = {lhs_val}")
    print(f"RHS = S({xi}) + S({eta}) = {s_xi} + {s_eta} = {rhs_val}\n")
    
    fkg_holds = lhs_val >= rhs_val
    print(f"Checking the condition: {lhs_val} >= {rhs_val} is {fkg_holds}.")
    
    if not fkg_holds:
        print("\nThe FKG lattice condition fails.")
        print("This counterexample on a graph with max degree d=1 shows that the positive")
        print("correlations property is not guaranteed for d>=1.")
        print("Thus, the largest integer d for which the property is guaranteed for ALL valid")
        print("graphs and ALL q>=2 must be 0.")

if __name__ == "__main__":
    main()