import math

def estimate_fatigue_life():
    """
    Estimates the fatigue life of a specimen using the Palmgren-Miner rule.
    """
    # Step 1: Define constants and problem parameters
    b = -0.09  # Basquin exponent
    
    # Assumed life in cycles at the endurance limit (sigma_e). This is a standard
    # engineering convention used when a specific value is not provided.
    N_e = 1_000_000
    
    # Stress levels as factors of the endurance limit sigma_e
    stress_factors = [1.0, 1.1, 1.2]
    
    # Fractions of total life spent at each stress level
    life_fractions = [0.70, 0.20, 0.10]
    
    # Step 2: Calculate cycles to failure (N_i) for each stress level
    # From Basquin's relation: N_i = N_e * (stress_factor_i)^(1/b)
    inv_b = 1.0 / b
    
    N_values = []
    for factor in stress_factors:
        N_i = N_e * (factor ** inv_b)
        N_values.append(N_i)
        
    N1, N2, N3 = N_values[0], N_values[1], N_values[2]

    # Step 3: Apply Palmgren-Miner rule
    # Sum(n_i / N_i) = 1, where n_i is the number of cycles applied at stress_i
    # n_i = life_fraction_i * N_total
    # N_total * Sum(life_fraction_i / N_i) = 1
    # N_total = 1 / Sum(life_fraction_i / N_i)
    
    damage_sum = 0
    for i in range(len(life_fractions)):
        damage_sum += life_fractions[i] / N_values[i]
        
    N_total = 1 / damage_sum

    # Step 4: Output the results as requested
    print("--- Fatigue Life Calculation ---")
    print(f"Basquin Exponent (b): {b}")
    print(f"Assumed Life at Endurance Limit (N_e): {int(N_e):,} cycles\n")

    print("Calculated cycles to failure for each stress level:")
    print(f"N1 (at 1.0 * sigma_e): {N1:,.0f} cycles")
    print(f"N2 (at 1.1 * sigma_e): {N2:,.0f} cycles")
    print(f"N3 (at 1.2 * sigma_e): {N3:,.0f} cycles\n")
    
    print("Final equation using Palmgren-Miner rule:")
    equation_str = f"N_total = 1 / ( {life_fractions[0]} / {N1:,.0f} + {life_fractions[1]} / {N2:,.0f} + {life_fractions[2]} / {N3:,.0f} )"
    print(equation_str)
    
    print("\n--- Final Answer ---")
    print(f"Estimated Total Fatigue Life (N_total): {N_total:.3e} cycles")


# Execute the function to get the answer
estimate_fatigue_life()
<<<4.710e+05>>>