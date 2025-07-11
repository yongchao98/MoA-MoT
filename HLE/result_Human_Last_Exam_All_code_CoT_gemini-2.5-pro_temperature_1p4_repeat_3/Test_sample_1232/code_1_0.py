import numpy as np
import itertools

def solve_potts_fkg():
    """
    This function numerically verifies the positive correlations property for the Potts model
    on a specific graph with maximum degree 3.

    The theoretical result is that the ferromagnetic Potts model exhibits the positive
    correlations property on ANY finite graph, regardless of its maximum degree.
    Therefore, there is no largest integer d for which this holds; it holds for all d.
    This corresponds to answer choice F.

    We demonstrate this principle by computing the covariance of two increasing functions
    for the Potts model on the complete graph K4, which has max_degree = 3.
    We expect the covariance to be non-negative.
    """

    # 1. Define graph and model parameters
    # Let's use the complete graph on 4 vertices, K_4.
    # The maximum degree of K_4 is 3.
    num_vertices = 4
    vertices = list(range(num_vertices))
    edges = list(itertools.combinations(vertices, 2))

    # Potts model parameters
    q = 3  # Number of states {1, 2, 3}
    beta = 0.5  # Inverse temperature (beta > 0 for ferromagnetic)

    # 2. Define two increasing functions
    # Let f be the spin at vertex 0 and g be the spin at vertex 1.
    # A function f(xi) = xi[v] is increasing because if xi <= eta, then by definition
    # xi[v] <= eta[v] for all v, which means f(xi) <= f(eta).
    v1 = 0
    v2 = 1

    print(f"Verifying the positive correlation property for the q={q} Potts model on K_4 (max_degree=3) with beta={beta}.")
    print(f"We will compute the covariance of the spins at vertex {v1} and vertex {v2}.")
    print(f"Let f(config) = config[{v1}] and g(config) = config[{v2}]. Both are increasing functions.")
    print("-" * 70)

    # 3. Iterate through all possible configurations to compute expectations
    partition_function = 0.0
    sum_f = 0.0
    sum_g = 0.0
    sum_fg = 0.0

    # Generate all possible configurations
    # A configuration is a mapping from V -> {1, ..., q}
    all_configs = itertools.product(range(1, q + 1), repeat=num_vertices)

    for config in all_configs:
        # Calculate the number of monochromatic edges, H(xi)
        num_monochromatic_edges = 0
        for u, v in edges:
            if config[u] == config[v]:
                num_monochromatic_edges += 1

        # Calculate the weight of the configuration
        weight = np.exp(2 * beta * num_monochromatic_edges)

        # Get the values of our functions for this configuration
        f_val = config[v1]
        g_val = config[v2]

        # Update the sums
        partition_function += weight
        sum_f += f_val * weight
        sum_g += g_val * weight
        sum_fg += f_val * g_val * weight

    # 4. Calculate expectations and covariance
    E_f = sum_f / partition_function
    E_g = sum_g / partition_function
    E_fg = sum_fg / partition_function
    covariance = E_fg - E_f * E_g

    # 5. Print the results of the equation Cov(f,g) = E[fg] - E[f] * E[g]
    print("The final equation is Cov(f,g) = E[f*g] - E[f]*E[g]\n")
    print(f"Calculated value for E[f*g]: {E_fg}")
    print(f"Calculated value for E[f]: {E_f}")
    print(f"Calculated value for E[g]: {E_g}")
    print(f"\nCov(f,g) = {E_fg} - {E_f} * {E_g}")
    print(f"Final Covariance: {covariance}")

    if covariance >= 0:
        print("\nAs expected, the covariance is non-negative.")
    else:
        print("\nWarning: The covariance is negative, which contradicts the established theorem.")
    
    print("-" * 70)
    print("This numerical check supports the theoretical result that the property holds for graphs with max_degree=3.")
    print("The general proof shows it holds for any max_degree.")
    print("Thus, there is no 'largest d'.")


solve_potts_fkg()
<<<F>>>