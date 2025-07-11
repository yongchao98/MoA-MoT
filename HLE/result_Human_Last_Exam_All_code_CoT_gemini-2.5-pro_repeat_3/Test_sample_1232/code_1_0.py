import numpy as np

def H_interaction(spin1, spin2):
    """
    Calculates the Potts model interaction term for a single edge.
    We ignore the constant factor 2*beta for simplicity, as it's positive
    and doesn't affect the supermodularity inequality.
    """
    return 1 if spin1 == spin2 else 0

def check_supermodularity(xi, eta):
    """
    Checks the supermodularity condition for two configurations on a two-vertex graph.
    xi = (spin_at_v1, spin_at_v2)
    eta = (spin_at_v1, spin_at_v2)
    """
    print(f"Let's test the supermodularity condition for the Potts model with q >= 3.")
    print(f"We check on a graph with two vertices and one edge (max degree = 1).")
    print(f"The condition is: H(xi v eta) + H(xi ^ eta) >= H(xi) + H(eta)")
    print("-" * 60)
    print(f"Configuration xi: {xi}")
    print(f"Configuration eta: {eta}")
    print("-" * 60)
    
    # Calculate the join (vee) and meet (wedge)
    xi_vee_eta = (max(xi[0], eta[0]), max(xi[1], eta[1]))
    xi_wedge_eta = (min(xi[0], eta[0]), min(xi[1], eta[1]))

    print(f"Component-wise maximum (vee): xi v eta = {xi_vee_eta}")
    print(f"Component-wise minimum (wedge): xi ^ eta = {xi_wedge_eta}")
    print("-" * 60)

    # Calculate the Hamiltonian for each configuration
    h_xi = H_interaction(xi[0], xi[1])
    h_eta = H_interaction(eta[0], eta[1])
    h_vee = H_interaction(xi_vee_eta[0], xi_vee_eta[1])
    h_wedge = H_interaction(xi_wedge_eta[0], xi_wedge_eta[1])

    # Check the inequality
    lhs = h_vee + h_wedge
    rhs = h_xi + h_eta

    print("Calculating each term in the inequality:")
    print(f"H(xi v eta) = I{{{xi_vee_eta[0]} = {xi_vee_eta[1]}}} = {h_vee}")
    print(f"H(xi ^ eta) = I{{{xi_wedge_eta[0]} = {xi_wedge_eta[1]}}} = {h_wedge}")
    print(f"H(xi)       = I{{{xi[0]} = {xi[1]}}} = {h_xi}")
    print(f"H(eta)      = I{{{eta[0]} = {eta[1]}}} = {h_eta}")
    print("-" * 60)
    
    print(f"Checking the inequality: {h_vee} + {h_wedge} >= {h_xi} + {h_eta}")
    print(f"Result: {lhs} >= {rhs}")

    if lhs >= rhs:
        print("\nThe supermodularity condition holds for this case.")
    else:
        print("\nThe supermodularity condition FAILS for this case.")

# Let's choose a counterexample for q=3 (spins {1, 2, 3})
# Let xi = (2, 2) and eta = (1, 3)
xi_config = (2, 2)
eta_config = (1, 3)
check_supermodularity(xi_config, eta_config)

print("\nConclusion:")
print("The supermodularity condition, which guarantees positive correlations, fails for a graph with maximum degree 1 when q >= 3.")
print("Since the property must hold for ANY q >= 2 and ANY graph with max_degree <= d, d cannot be 1 or greater.")
print("For d=0, the only connected graph is a single vertex with no edges. The property trivially holds.")
print("Therefore, the largest value for d is 0.")