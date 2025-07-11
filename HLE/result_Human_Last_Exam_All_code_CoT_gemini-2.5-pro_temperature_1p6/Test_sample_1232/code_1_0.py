import numpy as np

def check_fkg_for_potts():
    """
    This script demonstrates that the positive correlations property (FKG inequality)
    does not hold for the Potts model on a graph with max degree 1.
    We check the Holley criterion for a specific counterexample.
    """
    
    # 1. Define the counterexample parameters
    # Graph G: A single edge between two vertices. deg_max(G) = 1.
    # We only need to consider these two vertices.
    # Number of states q=3. Spins can be {1, 2, 3}.
    q = 3
    # Inverse temperature beta > 0. Let's use beta = 0.5 (so 2*beta = 1.0)
    beta = 0.5
    
    print("--- Counterexample Setup ---")
    print(f"Graph: A single edge (max degree d=1)")
    print(f"Number of states q: {q}")
    print(f"Inverse temperature beta: {beta}\n")

    # 2. Define the spin configurations xi and eta at the two vertices
    xi = np.array([2, 2])
    eta = np.array([1, 3])
    
    print("--- Spin Configurations ---")
    print(f"Configuration xi: {xi}")
    print(f"Configuration eta: {eta}\n")
    
    # 3. Calculate the combined configurations
    xi_max = np.maximum(xi, eta)
    xi_min = np.minimum(xi, eta)
    
    print("--- Combined Configurations ---")
    print(f"Configuration xi_vee_eta (pointwise max): {xi_max}")
    print(f"Configuration xi_wedge_eta (pointwise min): {xi_min}\n")

    # 4. Define the energy function H for a single edge
    # H = 2 * beta * I_{spin1 == spin2}
    def H(config):
        # I_{...} is 1 if spins are equal, 0 otherwise
        indicator = 1 if config[0] == config[1] else 0
        return 2 * beta * indicator

    # 5. Calculate the energy terms
    H_xi = H(xi)
    H_eta = H(eta)
    H_max = H(xi_max)
    H_min = H(xi_min)
    
    # The Holley criterion requires H(xi_max) + H(xi_min) >= H(xi) + H(eta)
    lhs = H_max + H_min
    rhs = H_xi + H_eta
    
    print("--- Checking the Holley Criterion ---")
    print("The FKG property holds if H(max) + H(min) >= H(xi) + H(eta)")
    
    print("\nCalculating individual H terms (H = 2*beta*I_{s1==s2}):")
    print(f"H(xi) for config {xi}: {H_xi:.1f} (since {xi[0]} == {xi[1]})")
    print(f"H(eta) for config {eta}: {H_eta:.1f} (since {eta[0]} != {eta[1]})")
    print(f"H(max) for config {xi_max}: {H_max:.1f} (since {xi_max[0]} != {xi_max[1]})")
    print(f"H(min) for config {xi_min}: {H_min:.1f} (since {xi_min[0]} != {xi_min[1]})")

    print("\nFinal Inequality Check:")
    print(f"  H(max) + H(min)  >=  H(xi) + H(eta)")
    # The final print statement requested by the prompt
    print(f"   {H_max:.1f}   +   {H_min:.1f}   >=   {H_xi:.1f}  +   {H_eta:.1f}")
    print(f"        {lhs:.1f}          >=         {rhs:.1f}")
    print(f"\nThis inequality is {lhs >= rhs}.")
    
    print("\n--- Conclusion ---")
    print("Since the property fails for a graph with max degree d=1,")
    print("the largest integer d for which it holds for *all* cases must be d=0.")

check_fkg_for_potts()