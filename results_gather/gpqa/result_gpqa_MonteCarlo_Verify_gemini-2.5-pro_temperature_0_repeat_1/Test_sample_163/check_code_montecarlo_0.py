import numpy as np

def solve_binary_star_mass_ratio():
    """
    Solves the binary star mass ratio problem using both Monte Carlo
    simulation and deterministic calculation.
    """
    # --- (a) Sampling (Monte Carlo Exploration) ---
    print("--- Monte Carlo Exploration ---")
    
    # Given mean values
    params = {
        'P1': 2.0,    # years
        'K1_1': 10.0, # km/s
        'K1_2': 5.0,  # km/s
        'P2': 1.0,    # years
        'K2_1': 15.0, # km/s
        'K2_2': 10.0  # km/s
    }
    
    num_samples = 100000
    # Assume a small relative uncertainty (e.g., 0.5%) for simulation
    relative_uncertainty = 0.005
    
    # Generate samples from normal distributions
    samples = {key: np.random.normal(loc=val, scale=val*relative_uncertainty, size=num_samples)
               for key, val in params.items()}

    # Calculate the mass ratio for each sample
    # M1/M2 = (P1/P2) * ((K1_1 + K1_2) / (K2_1 + K2_2))^3
    k_sum_1 = samples['K1_1'] + samples['K1_2']
    k_sum_2 = samples['K2_1'] + samples['K2_2']
    
    sampled_ratios = (samples['P1'] / samples['P2']) * (k_sum_1 / k_sum_2)**3
    
    # --- (b) Narrowing Candidates ---
    mean_ratio = np.mean(sampled_ratios)
    std_ratio = np.std(sampled_ratios)
    
    print(f"Simulation with {num_samples} samples suggests a mass ratio.")
    print(f"Mean Ratio: {mean_ratio:.4f}")
    print(f"Standard Deviation: {std_ratio:.4f}")
    print(f"This strongly suggests the answer is around {mean_ratio:.2f}.")
    
    # --- (c) Exact Checks (Deterministic Verification) ---
    print("\n--- Deterministic Verification ---")
    
    P1 = params['P1']
    P2 = params['P2']
    K1_sum = params['K1_1'] + params['K1_2'] # 10 + 5 = 15 km/s
    K2_sum = params['K2_1'] + params['K2_2'] # 15 + 10 = 25 km/s
    
    # The formula for the ratio of total masses M1/M2
    exact_ratio = (P1 / P2) * (K1_sum / K2_sum)**3
    
    print(f"System 1 Period (P1): {P1} years")
    print(f"System 2 Period (P2): {P2} year")
    print(f"System 1 Sum of RV Amplitudes (K1_sum): {K1_sum} km/s")
    print(f"System 2 Sum of RV Amplitudes (K2_sum): {K2_sum} km/s")
    print(f"Ratio calculation: ({P1}/{P2}) * ({K1_sum}/{K2_sum})**3")
    print(f"Exact Ratio M1/M2 = {exact_ratio:.4f}")
    
    # --- (d) Print Final Verified Answer ---
    options = {'A': 0.7, 'B': 0.4, 'C': 0.6, 'D': 1.2}
    
    # Find the closest option to the exact result
    best_option = min(options, key=lambda k: abs(options[k] - exact_ratio))
    
    print("\n--- Final Answer ---")
    print(f"The calculated ratio is {exact_ratio:.4f}, which is closest to option {best_option}'s value of {options[best_option]}.")
    print(f"Therefore, system_1 is approximately {exact_ratio:.3f} times the mass of system_2.")
    print(f"The correct option is {best_option}.")

# Run the solver
solve_binary_star_mass_ratio()