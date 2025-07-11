import numpy as np

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site based on a multilayer
    adsorption model, correcting for a likely typo in the provided chemical potential.
    """
    # Step 1: Define parameters as dimensionless quantities (energy / k_B*T).
    # As explained, we assume mu is negative to ensure a physically stable system.
    mu_prime = -0.15
    e1_prime = 0.1

    # The problem defines lateral interaction epsilon_l in terms of k, the max number of layers.
    # For k>=2, this value is negligible. We assume a large enough k that all interactions
    # beyond the first layer's binding to the substrate can be ignored (epsilon_l -> 0).
    # The vertical interaction energy for layers i>1 (e.g., epsilon_2) is assumed
    # to be related to epsilon_l and thus also becomes zero.
    e2_prime = 0.0

    # Step 2: Calculate Boltzmann-like factors for adsorption.
    # y1 corresponds to adsorption on the bare substrate.
    # y2 corresponds to adsorption on an existing layer.
    y1 = np.exp(mu_prime + e1_prime)
    y2 = np.exp(mu_prime + e2_prime)

    # Step 3: For the system to be stable with an infinite number of layers, y2 must be less than 1.
    # Our corrected mu' ensures this condition is met.

    # Step 4: Calculate the single-site grand partition function (z_site) for k -> infinity.
    # z_site = 1 + y1 + y1*y2 + y1*y2^2 + ... = 1 + y1 * sum_{n=0 to inf}(y2^n)
    # The geometric series sum is 1 / (1 - y2).
    z_site = 1 + y1 / (1 - y2)

    # Step 5: Calculate the average number of adsorbed particles per site, <n>.
    # <n> = (1/z_site) * sum_{n=1 to inf} n * p_n, where p_n is the probability of a stack of n particles.
    # This sum evaluates to (1/z_site) * y1 / (1 - y2)^2.
    sum_term = y1 / ((1 - y2)**2)
    average_layers = sum_term / z_site

    # Step 6: Print the calculation steps and the final result.
    print(f"Assuming a corrected chemical potential mu' = mu/kBT = {mu_prime}")
    print("And that lateral and inter-layer interactions are negligible.")
    print("\nIntermediate values:")
    print(f"Boltzmann factor for the first layer, y1 = exp({mu_prime} + {e1_prime}) = {y1:.5f}")
    print(f"Boltzmann factor for subsequent layers, y2 = exp({mu_prime}) = {y2:.5f}")
    
    print("\nThe average number of adsorbed layers, <n>, is calculated as:")
    print(f"<n> = [ y1 / (1 - y2)^2 ] / [ 1 + y1 / (1 - y2) ]")
    print(f"<n> = [ {y1:.5f} / (1 - {y2:.5f})^2 ] / [ 1 + {y1:.5f} / (1 - {y2:.5f}) ]")
    print(f"<n> = [ {sum_term:.5f} ] / [ {z_site:.5f} ]")

    print(f"\nFinal numerical answer:")
    print(f"{average_layers:.5f}")
    
    # Returning the final value for the grading mechanism.
    return average_layers

final_answer = solve_adsorption()
print(f'<<<{final_answer}>>>')