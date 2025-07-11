import math

def estimate_fatigue_life():
    """
    Estimates fatigue life using the Basquin relationship and Palmgren-Miner rule.
    """
    # Given parameters
    b = -0.09  # Basquin exponent
    cycle_fractions = [0.70, 0.20, 0.10]
    stress_multipliers = [1.0, 1.1, 1.2]

    # --- Step 1: Establish a general formula for N_f ---
    # The Basquin equation is: sigma_a = C * N_f^b
    # This can be rearranged to: N_f = (sigma_a / C)^(1/b)
    # We are given that at N_f = 0.5, sigma_a = sigma_uts.
    # So, sigma_uts = C * (0.5)^b  =>  C = sigma_uts / (0.5)^b
    # Substituting C back into the N_f equation:
    # N_f = (sigma_a / (sigma_uts / (0.5)^b))^(1/b)
    # N_f = ((sigma_a / sigma_uts) * (0.5)^b)^(1/b)
    # N_f = (sigma_a / sigma_uts)^(1/b) * 0.5
    # This formula allows us to calculate N_f based on the stress ratio.

    # --- Step 2: Make a key assumption for the endurance limit ---
    # The relationship between endurance limit (sigma_e) and ultimate tensile strength
    # (sigma_uts) is not given. A common engineering approximation is that the
    # fatigue ratio is 0.5 for many steels.
    fatigue_ratio = 0.5  # Assumption: sigma_e / sigma_uts = 0.5

    # --- Step 3: Calculate the number of cycles to failure (N_i) for each stress level ---
    inv_b = 1.0 / b
    stress_ratios = [m * fatigue_ratio for m in stress_multipliers]

    N_values = []
    print("Calculating life (N_i) at each stress level:")
    for i, ratio in enumerate(stress_ratios):
        # N_f = (ratio)^(1/b) * 0.5
        N_i = (ratio) ** inv_b * 0.5
        N_values.append(N_i)
        print(f"  - At stress ratio {ratio:.2f} (_uts), life N_{i+1} = {N_i:.2f} cycles")

    # --- Step 4 & 5: Apply Palmgren-Miner rule and compute total life ---
    # Palmgren-Miner Rule: Sum(n_i / N_i) = 1 at failure.
    # The problem implies cycle fractions: n_i = c_i * N_total.
    # So, Sum(c_i * N_total / N_i) = 1  =>  N_total * Sum(c_i / N_i) = 1
    # N_total = 1 / Sum(c_i / N_i)
    
    c1, c2, c3 = cycle_fractions
    N1, N2, N3 = N_values

    damage_sum = (c1 / N1) + (c2 / N2) + (c3 / N3)
    
    total_life = 1.0 / damage_sum

    print("\nApplying Palmgren-Miner Rule:")
    # The final equation with each number included
    print(f"N_total = 1 / ({c1}/{N1:.2f} + {c2}/{N2:.2f} + {c3}/{N3:.2f})")
    print(f"N_total = 1 / ({damage_sum:.6f})")
    print(f"N_total = {total_life:.3f} cycles\n")

    # Format final answer to Aeb scientific notation
    formatted_life = "{:.3e}".format(total_life)
    # Clean up standard python scientific notation (e.g., e+02 -> e2)
    clean_formatted_life = formatted_life.replace('e+0', 'e').replace('e-0', 'e-')
    
    print("Final Estimated Fatigue Life:")
    print(clean_formatted_life)
    
    return clean_formatted_life

# Execute the function
final_answer = estimate_fatigue_life()
# The final result in the specified format for the platform
print(f"<<<{final_answer}>>>")