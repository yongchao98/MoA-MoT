import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site based on a
    finite-layer adsorption model derived from the problem description.
    """
    # Parameters derived from the problem statement and our assumptions.
    # Energies are in units of k_B*T, so we can treat k_B*T as 1.
    mu = 0.15
    eps1 = 0.1
    eps2 = 0.02  # Assumed from a typo in ε_l
    k_max = 4     # Inferred from z_l and z_inter

    # Calculate the statistical weights P(L) for a site having L layers.
    # P(L) is proportional to exp(-β(E_L - L*μ)).
    # We set P(0) = 1 as the reference weight.
    P = [0.0] * (k_max + 1)
    P[0] = 1.0

    # C = exp(μ - ε_1) corresponds to the weight of adding the first layer.
    C = math.exp(mu - eps1)
    # X = exp(μ - ε_2) corresponds to the weight factor for each subsequent layer.
    X = math.exp(mu - eps2)

    # For L > 0, P(L) = C * X^(L-1)
    for L in range(1, k_max + 1):
        P[L] = C * (X**(L - 1))

    # The numerator for the average is the sum of L * P(L)
    numerator_terms = [L * P[L] for L in range(k_max + 1)]
    numerator = sum(numerator_terms)

    # The denominator is the single-site partition function, which is the sum of all weights.
    denominator = sum(P)

    # Calculate the average number of layers
    avg_layers = numerator / denominator

    # Print the final equation with all its numerical components
    print("The average number of layers <k> is calculated using the formula:")
    print(f"<k> = Σ(L * P(L)) / Σ(P(L)) for L = 0 to {k_max}\n")

    print("The individual statistical weights P(L) are:")
    for L in range(k_max + 1):
        print(f"P({L}) = {P[L]:.5f}")

    print("\nBreaking down the formula with these values:")
    numerator_str = " + ".join([f"{val:.5f}" for val in numerator_terms])
    denominator_str = " + ".join([f"{p:.5f}" for p in P])
    
    print(f"<k> = ({numerator_str}) / ({denominator_str})")
    print(f"<k> = {numerator:.5f} / {denominator:.5f}")
    
    print(f"\nThe calculated average number of adsorbed layers per site is: {avg_layers:.4f}")

    # Return the final answer in the specified format
    print(f"\n<<<{avg_layers:.4f}>>>")

solve_adsorption()