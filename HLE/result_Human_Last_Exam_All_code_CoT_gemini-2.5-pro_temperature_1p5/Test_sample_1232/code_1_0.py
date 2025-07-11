import numpy as np

def calculate_hamiltonian(config, beta):
    """
    Calculates the Potts Hamiltonian for a 2-vertex, 1-edge graph.
    H(config) = 2 * beta * I_{config[0] == config[1]}
    """
    v1_spin, v2_spin = config
    # Indicator function I_{spin1 == spin2}
    indicator = 1 if v1_spin == v2_spin else 0
    return 2 * beta * indicator

def main():
    """
    Checks the Holley submodularity condition for a specific counterexample.
    """
    # We can choose any beta > 0. Let's use beta = 1.0
    beta = 1.0
    # This counterexample works for any q >= 3.
    q = 3

    print(f"Checking Holley condition for the Potts model with q={q}, beta={beta}")
    print("Graph: Two vertices connected by one edge (max_degree = 1).\n")

    # Define the two configurations xi and eta
    xi = (2, 2)
    eta = (1, 3)

    print(f"Configuration xi: {xi}")
    print(f"Configuration eta: {eta}\n")

    # Define the meet and join configurations
    xi_meet_eta = (min(xi[0], eta[0]), min(xi[1], eta[1]))
    xi_join_eta = (max(xi[0], eta[0]), max(xi[1], eta[1]))

    print(f"Meet (xi_wedge_eta): {xi_meet_eta}")
    print(f"Join (xi_vee_eta):   {xi_join_eta}\n")

    # Calculate the Hamiltonian for each configuration
    H_xi = calculate_hamiltonian(xi, beta)
    H_eta = calculate_hamiltonian(eta, beta)
    H_meet = calculate_hamiltonian(xi_meet_eta, beta)
    H_join = calculate_hamiltonian(xi_join_eta, beta)
    
    # Check the inequality: H(meet) + H(join) >= H(xi) + H(eta)
    lhs = H_meet + H_join
    rhs = H_xi + H_eta
    is_satisfied = lhs >= rhs

    print("The submodularity condition is: H(meet) + H(join) >= H(xi) + H(eta)")
    print("Let's plug in the numbers:")
    print(f"Equation: {H_meet:.1f} + {H_join:.1f} >= {H_xi:.1f} + {H_eta:.1f}")
    print(f"Result:   {lhs:.1f} >= {rhs:.1f}")
    print(f"\nIs the condition satisfied? {is_satisfied}")
    
    if not is_satisfied:
        print("\nSince the condition fails for a graph with max_degree=1 (for q>=3),")
        print("the positive correlations property does not generally hold for d=1.")
        print("This implies it also doesn't hold for any d > 1.")

    print("\nFor d=0, the graph is a single vertex, so there are no edges.")
    print("The Hamiltonian is always 0, and the measure is uniform.")
    print("The positive correlation property holds in this case due to Chebyshev's sum inequality.")
    print("\nTherefore, the largest integer d for which the property always holds is 0.")


if __name__ == "__main__":
    main()
