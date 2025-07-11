import numpy as np

def potts_hamiltonian(config, edges, beta):
    """Calculates the Hamiltonian for a given configuration."""
    energy = 0
    for u, v in edges:
        if config[u] == config[v]:
            energy += 1
    return 2 * beta * energy

# Define the parameters for the counterexample
# 1. Graph G = P_2 (two vertices, one edge). deg_max(G) = 1.
V = [0, 1]
edges = [(0, 1)]

# 2. Number of states q = 3.
# 3. Inverse temperature beta > 0.
beta = 1.0

# 4. Two configurations xi and eta.
# We use 1-based indexing for states {1, 2, 3} as in the problem.
xi = np.array([2, 2])
eta = np.array([1, 3])

# Calculate the meet (min) and join (max) configurations
xi_min = np.minimum(xi, eta)
xi_max = np.maximum(xi, eta)

print(f"Demonstrating a counterexample for d=1, q=3, beta={beta}")
print("-" * 60)
print(f"Configuration xi: {xi}")
print(f"Configuration eta: {eta}")
print(f"Meet (xi_min): {xi_min}")
print(f"Join (xi_max): {xi_max}")
print("-" * 60)

# Calculate the Hamiltonian for each configuration
H_xi = potts_hamiltonian(xi, edges, beta)
H_eta = potts_hamiltonian(eta, edges, beta)
H_xi_min = potts_hamiltonian(xi_min, edges, beta)
H_xi_max = potts_hamiltonian(xi_max, edges, beta)

# Check Holley's criterion: H(xi_min) + H(xi_max) >= H(xi) + H(eta)
lhs = H_xi_min + H_xi_max
rhs = H_xi + H_eta

print("Checking Holley's criterion: H(xi_min) + H(xi_max) >= H(xi) + H(eta)")
print(f"H(xi) = 2 * {beta} * I({xi[0]} == {xi[1]}) = {H_xi}")
print(f"H(eta) = 2 * {beta} * I({eta[0]} == {eta[1]}) = {H_eta}")
print(f"H(xi_min) = 2 * {beta} * I({xi_min[0]} == {xi_min[1]}) = {H_xi_min}")
print(f"H(xi_max) = 2 * {beta} * I({xi_max[0]} == {xi_max[1]}) = {H_xi_max}")
print("-" * 60)

# Print the final equation with numbers
print("Final inequality check:")
print(f"LHS = H(xi_min) + H(xi_max) = {H_xi_min} + {H_xi_max} = {lhs}")
print(f"RHS = H(xi) + H(eta) = {H_xi} + {H_eta} = {rhs}")
print(f"Is {lhs} >= {rhs}? {lhs >= rhs}")
print("\nSince the inequality is false, the positive correlations property does not hold for this case.")
print("This proves the statement is false for d=1.")
