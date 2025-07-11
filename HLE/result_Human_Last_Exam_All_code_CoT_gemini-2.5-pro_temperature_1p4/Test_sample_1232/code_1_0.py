import numpy as np
import itertools

def solve_potts_correlation():
    """
    Numerically investigates the positive correlation property for the Potts model
    on a specific graph.
    """
    # Parameters for the Potts model
    q = 3  # Number of states. The FKG lattice condition fails for q>=3.
    beta = 1.0  # Inverse temperature
    n_vertices = 3 # Number of vertices in the graph
    
    # Define the graph as a list of edges. This is a path graph P_3: 0--1--2.
    # The maximum degree of this graph is d=2.
    edges = [(0, 1), (1, 2)]

    # Variables to store the sums needed for expectations
    Z = 0.0      # Normalization constant (partition function)
    sum_f = 0.0    # Sum of f(xi) * weight(xi)
    sum_g = 0.0    # Sum of g(xi) * weight(xi)
    sum_fg = 0.0   # Sum of f(xi) * g(xi) * weight(xi)

    # The spin states are {1, 2, ..., q}
    states = range(1, q + 1)

    # Iterate over all possible configurations xi in {1,...,q}^V
    for xi in itertools.product(states, repeat=n_vertices):
        # A configuration xi is a tuple, e.g., (1, 2, 1)
        
        # Calculate the energy term H(xi) = sum I(xi(x) = xi(y))
        hamiltonian = 0
        for u, v in edges:
            if xi[u] == xi[v]:
                hamiltonian += 1
        
        # Calculate the Gibbs weight for the configuration
        weight = np.exp(2 * beta * hamiltonian)
        
        # Define two increasing functions, f and g.
        # f(xi) = spin of the first vertex
        # g(xi) = spin of the third vertex
        # These are increasing because a larger spin value at a vertex
        # leads to a larger (or equal) configuration in the coordinatewise order.
        f_val = xi[0]
        g_val = xi[2]
        
        # Accumulate the weighted sums
        Z += weight
        sum_f += f_val * weight
        sum_g += g_val * weight
        sum_fg += f_val * g_val * weight

    # Calculate the final expectations by normalizing
    E_f = sum_f / Z
    E_g = sum_g / Z
    E_fg = sum_fg / Z

    # Calculate the covariance
    covariance = E_fg - E_f * E_g

    # Print the results as an equation
    print("Checking the positive correlation property: E[fg] >= E[f]E[g]")
    print(f"Graph: Path on {n_vertices} vertices (max degree = {len(edges)})")
    print(f"Parameters: q = {q}, beta = {beta}")
    print(f"We test with increasing functions f(xi) = xi[0] and g(xi) = xi[2].")
    print("\nCalculating the final equation: E[fg] - E[f] * E[g]")
    print(f"{E_fg:.6f} - {E_f:.6f} * {E_g:.6f} = {covariance:.6f}")
    
    if covariance >= 0:
        print("\nThe result is non-negative, as predicted by the theory.")
        print("This supports the conclusion that the property holds regardless of max degree.")
    else:
        print("\nA counter-example was found!")

solve_potts_correlation()
<<<F>>>