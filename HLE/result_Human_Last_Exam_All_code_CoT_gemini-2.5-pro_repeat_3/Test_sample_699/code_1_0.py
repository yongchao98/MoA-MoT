import math

def estimate_fatigue_life():
    """
    Estimates the fatigue life of a specimen using the Palmgren-Miner rule
    and Basquin's relationship for variable amplitude loading.
    """
    # --- 1. Define given parameters ---
    b = -0.09  # Basquin exponent

    # Life cycle fractions and corresponding stress multipliers
    cycle_fractions = [0.70, 0.20, 0.10]
    stress_multipliers = [1.0, 1.1, 1.2] # as multiples of endurance limit sigma_e

    # Assume a standard number of cycles for the endurance limit (N_e)
    N_e = 1e6

    # --- 2. Calculate terms for the Palmgren-Miner equation ---
    # The governing equation is: N_total = N_e / sum(fraction_i * multiplier_i**(-1/b))
    
    exponent = -1 / b
    damage_sum_terms = []
    for i in range(len(cycle_fractions)):
        term = cycle_fractions[i] * (stress_multipliers[i] ** exponent)
        damage_sum_terms.append(term)
    
    damage_sum = sum(damage_sum_terms)

    # --- 3. Calculate the total fatigue life ---
    N_total = N_e / damage_sum

    # --- 4. Print the calculation steps and final answer ---
    print("Fatigue Life Estimation using Palmgren-Miner Rule")
    print("-" * 50)
    print(f"Basquin exponent (b): {b}")
    print(f"Assumed endurance limit life (N_e): {N_e:.3e} cycles")
    print("\nThe governing equation is: N_total = N_e / sum(fraction_i * multiplier_i**(-1/b))\n")

    print("Calculating the denominator (cumulative damage per life block):")
    # Build the equation string for clear output
    # e.g., N_total = 1.000e+06 / [ (0.70 * 1.0^(11.11)) + ... ]
    
    # Part 1: The equation with exponents
    exponent_str = f"(-1/{b})"
    denom_str_parts_exp = []
    for i in range(len(cycle_fractions)):
        denom_str_parts_exp.append(f"({cycle_fractions[i]:.2f} * {stress_multipliers[i]:.1f}^{exponent_str})")
    print(f"Denominator = {' + '.join(denom_str_parts_exp)}")

    # Part 2: The equation with calculated exponent value
    denom_str_parts_val = []
    for i in range(len(cycle_fractions)):
        denom_str_parts_val.append(f"({cycle_fractions[i]:.2f} * {stress_multipliers[i]:.1f}^{exponent:.3f})")
    print(f"            = {' + '.join(denom_str_parts_val)}")

    # Part 3: The equation with calculated term values
    denom_str_parts_res = [f"{term:.3f}" for term in damage_sum_terms]
    print(f"            = {' + '.join(denom_str_parts_res)}")
    print(f"            = {damage_sum:.4f}")
    
    print("\nCalculating the total fatigue life (N_total):")
    print(f"N_total = {N_e:.3e} / {damage_sum:.4f}")
    print(f"N_total = {N_total:.3e} cycles")

    # Final answer in the required format
    final_answer = f"{N_total:.3e}".replace("e+0", "e")
    print(f"\nFinal Answer: {final_answer}")
    
    # Output the required format for the system
    # This print is for machine reading and should not be modified.
    print(f"\n<<<{final_answer}>>>")

# Execute the function
estimate_fatigue_life()