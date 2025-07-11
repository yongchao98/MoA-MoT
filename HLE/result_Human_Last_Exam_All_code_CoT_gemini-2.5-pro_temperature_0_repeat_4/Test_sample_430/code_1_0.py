import numpy as np

def solve_adsorption_problem():
    """
    Solves the multi-layer adsorption problem based on a set of assumptions
    to resolve ambiguities in the problem statement.
    """
    # Parameters from the problem statement, in units of k_B*T
    # beta_e1 is the dimensionless energy for the first layer (epsilon_1 / k_B*T)
    beta_e1 = 0.1
    # beta_mu is the dimensionless chemical potential (mu / k_B*T)
    beta_mu = 0.15
    # z_l is the lateral coordination number
    zl = 4
    # z_inter is the vertical coordination number
    z_inter = 4
    # T is the temperature in Kelvin
    T = 318

    # --- Assumptions to make the problem solvable ---
    # 1. The maximum number of layers, k_max, is assumed to be equal to z_inter.
    k_max = z_inter

    # 2. The lateral interaction energy is given by beta_el = (0.02)^k.
    # For k_max = 4, beta_el = (0.02)**4 = 1.6e-7, which is negligible.
    # We assume el = 0, which simplifies the mean-field calculation.
    beta_el = 0.0

    # 3. The interaction energy for subsequent layers (j > 1), beta_ev, is not given.
    # We assume it is equal to the chemical potential, beta_mu, as this uses a
    # given parameter to resolve the ambiguity.
    beta_ev = beta_mu

    print("--- Model Parameters and Assumptions ---")
    print(f"Max layers (k_max) = {k_max}")
    print(f"beta*epsilon_1 = {beta_e1}")
    print(f"beta*mu = {beta_mu}")
    print(f"beta*epsilon_j (for j>1) = {beta_ev} (assumed)")
    print(f"beta*epsilon_l = {beta_el} (negligible)\n")


    # --- Calculation ---
    # The model simplifies to independent sites with no lateral interactions.
    # The fugacity-like factor for adding a particle to layer j is:
    # X_j = exp(beta * (mu + ej))

    X1 = np.exp(beta_mu + beta_e1)
    Xv = np.exp(beta_mu + beta_ev)  # For j > 1

    # Calculate the unnormalized probability (term in z_site) for having m layers.
    # Term_0 corresponds to an empty site (energy=0), so its Boltzmann factor is 1.
    # Term_m = product of X_j from j=1 to m.
    terms = [1.0]
    current_term = 1.0
    for m in range(1, k_max + 1):
        if m == 1:
            current_term *= X1
        else:
            current_term *= Xv
        terms.append(current_term)

    # The single-site partition function, z_site, is the sum of these terms.
    z_site = sum(terms)

    # The probability of having exactly m layers is P_m = Term_m / z_site.
    probabilities = [term / z_site for term in terms]

    # The average number of layers, <k>, is the expectation value of m: Sum(m * P_m).
    n_avg = sum(m * p for m, p in enumerate(probabilities))

    # --- Output ---
    print("--- Results ---")
    print(f"Single-site partition function z = {z_site:.4f}")
    print("Probabilities for m = 0, 1, 2, 3, 4 layers:")
    for m, p in enumerate(probabilities):
        print(f"P(m={m}) = {p:.4f}")

    print("\nThe final equation for the average number of layers <k> is:")
    equation_parts = [f"({m} * {probabilities[m]:.4f})" for m in range(1, k_max + 1)]
    equation_str = "<k> = " + " + ".join(equation_parts)
    print(equation_str)

    print(f"\nThe calculated average number of adsorbed layers per site is:")
    print(f"<k> = {n_avg:.4f}")

solve_adsorption_problem()