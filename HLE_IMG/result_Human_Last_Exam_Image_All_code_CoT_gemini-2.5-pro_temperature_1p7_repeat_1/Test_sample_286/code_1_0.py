import numpy as np

def calculate_entropy(r):
    """Calculates the von Neumann entropy for a given Bloch vector length r."""
    # Ensure r is within the valid range [0, 1]
    r = np.clip(r, 0.0, 1.0)
    # Avoid log(0) for pure states
    if np.isclose(r, 1.0):
        return 0.0
    
    # Eigenvalues of the density matrix
    p1 = (1 + r) / 2
    p2 = (1 - r) / 2
    
    # Von Neumann entropy formula
    entropy = -(p1 * np.log2(p1) + p2 * np.log2(p2))
    return entropy

# --- Analysis of Plot B ---

# Values read from the graph at time t=0
sz_0 = 0.5
splus_0 = 0.7
S_graph_0 = 0.0

# Values read from the graph at time t=10
sz_10 = 0.6
splus_10 = 0.7
S_graph_10 = 0.25

# Calculate the squared Bloch vector length r^2
r_sq_0 = sz_0**2 + splus_0**2
r_sq_10 = sz_10**2 + splus_10**2

print("Analysis of Plot B shows a physical contradiction:")
print(f"At t=0: <σz>≈{sz_0}, |<σ+>|≈{splus_0}. This gives r^2 ≈ {r_sq_0:.2f}.")
print(f"At t=10: <σz>≈{sz_10}, |<σ+>|≈{splus_10}. This gives r^2 ≈ {r_sq_10:.2f}.")
print("")
print(f"The graph shows entropy S increasing from {S_graph_0} to {S_graph_10}.")
print("An increase in entropy signifies that the state is becoming more mixed.")
print("Therefore, the length of the Bloch vector squared (r^2) must DECREASE.")
print(f"However, the values from the plot show that r^2 INCREASED from {r_sq_0:.2f} to {r_sq_10:.2f}.")
print("This is a fundamental contradiction, so Plot B cannot be a valid quantum evolution.")
