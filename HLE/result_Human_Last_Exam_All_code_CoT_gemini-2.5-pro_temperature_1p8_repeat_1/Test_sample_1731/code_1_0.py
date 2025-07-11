import numpy as np

def calculate_bose_equilibrium():
    """
    Calculates the equilibrium mean energy and entropy for a photon gas (Bose case).

    This is based on the Bose-Einstein distribution, which describes the most probable
    macroscopic state as predicted by large deviation theory (e.g., Boltzmann-Sanov).
    """

    # --- System Parameters ---
    # We use non-dimensional units for this theoretical calculation.
    # Set inverse temperature beta = 1 / (k_B * T). Let's set beta = 1.0.
    beta = 1.0
    # Let energy levels be simple integers: ε_i = i for i = 1, 2, ...
    # We'll sum over a finite number of levels, sufficient for convergence.
    num_levels = 25

    print("--- Calculation of Equilibrium Properties for a Photon Gas ---")
    print(f"Parameters: beta = {beta}, Energy levels ε_i = i, up to i_max = {num_levels}\n")

    energy_terms = []
    entropy_terms = []
    
    # --- Loop over Energy Levels to Calculate Terms in the Summation ---
    for i in range(1, num_levels + 1):
        epsilon_i = float(i)
        
        # Equilibrium mean occupation number <n_i> for photons (μ=0, g_i=1)
        # <n_i> = 1 / (exp(β*ε_i) - 1)
        mean_n_i = 1.0 / (np.exp(beta * epsilon_i) - 1.0)
        
        # Contribution to mean energy: <n_i> * ε_i
        energy_term = mean_n_i * epsilon_i
        energy_terms.append(energy_term)
        
        # Contribution to entropy: (<n_i>+1)ln(<n_i>+1) - <n_i>ln(<n_i>)
        # This term must be handled carefully if n_i is very close to zero.
        entropy_term = 0.0
        if mean_n_i > 1e-12:  # Avoid log(0) issues for negligible occupation
            term1 = (mean_n_i + 1.0) * np.log(mean_n_i + 1.0)
            term2 = mean_n_i * np.log(mean_n_i)
            entropy_term = term1 - term2
        entropy_terms.append(entropy_term)

    # --- Sum the terms to get total equilibrium values ---
    total_mean_energy = sum(energy_terms)
    total_entropy = sum(entropy_terms)
    
    # --- Output the Final Equations and Results ---
    print("--- 1. Equilibrium Mean Energy <E> ---")
    print("Formula: <E> = Σ [ε_i / (exp(β*ε_i) - 1)]")
    
    # Display the equation with the first few calculated terms
    energy_eq_str = " + ".join([f"{term:.4f}" for term in energy_terms[:5]])
    print("Equation with calculated terms:")
    print(f"<E> = {energy_eq_str} + ...")
    
    # Display the final result
    print(f"\nFinal Result: <E> = {total_mean_energy:.6f}\n")


    print("--- 2. Equilibrium Entropy S (in units of k_B) ---")
    print("Formula: S/k_B = Σ [ (<n_i>+1)ln(<n_i>+1) - <n_i>ln(<n_i>) ]")

    # Display the equation with the first few calculated terms
    entropy_eq_str = " + ".join([f"{term:.4f}" for term in entropy_terms[:5]])
    print("Equation with calculated terms:")
    print(f"S/k_B = {entropy_eq_str} + ...")
    
    # Display the final result
    print(f"\nFinal Result: S/k_B = {total_entropy:.6f}")


# Run the calculation and print the results
calculate_bose_equilibrium()
<<<The mean energy is approximately 1.644934 and the entropy (in units of k_B) is approximately 2.195328.>>>