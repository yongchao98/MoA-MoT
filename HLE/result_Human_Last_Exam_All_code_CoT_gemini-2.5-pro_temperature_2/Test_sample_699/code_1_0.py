import math

def estimate_fatigue_life():
    """
    Estimates the fatigue life of a specimen using Basquin's relationship and the Palmgren-Miner rule.
    """
    # Given parameters
    b = -0.09  # Basquin exponent

    # Loading history fractions
    frac1 = 0.70
    frac2 = 0.20
    frac3 = 0.10

    # Step 1 & 2: Define stress levels as a ratio of ultimate tensile strength (sigma_u).
    # We assume the endurance limit (sigma_e) is 0.5 * sigma_u.
    stress_ratio_1 = 0.50  # for sigma_e
    stress_ratio_2 = 1.1 * stress_ratio_1  # for 1.1 * sigma_e
    stress_ratio_3 = 1.2 * stress_ratio_1  # for 1.2 * sigma_e

    # Step 3: Calculate the number of cycles to failure (N_i) for each stress level.
    # The formula is N_i = (stress_ratio_i)^(1/b) * 0.5
    # This is derived from N = (sigma_a / C)^(1/b) where C = sigma_u / (0.5)^b
    try:
        inv_b = 1.0 / b
        N1 = math.pow(stress_ratio_1, inv_b) * 0.5
        N2 = math.pow(stress_ratio_2, inv_b) * 0.5
        N3 = math.pow(stress_ratio_3, inv_b) * 0.5
    except (ValueError, ZeroDivisionError) as e:
        print(f"An error occurred during calculation: {e}")
        return

    # Step 4: Apply the Palmgren-Miner linear damage rule.
    # Total Damage D = sum(n_i / N_i) = 1.
    # N_f * (frac1/N1 + frac2/N2 + frac3/N3) = 1
    # N_f = 1 / (frac1/N1 + frac2/N2 + frac3/N3)
    damage_sum = (frac1 / N1) + (frac2 / N2) + (frac3 / N3)
    if damage_sum == 0:
        print("Cannot calculate fatigue life, total damage is zero.")
        return
        
    total_fatigue_life = 1.0 / damage_sum

    # Print the details of the calculation
    print("--- Fatigue Life Estimation ---")
    print(f"Basquin Exponent (b): {b}\n")
    print("Cycles to failure for each stress level:")
    print(f"N1 (at {stress_ratio_1:.2f}*sigma_u): {N1:.3f} cycles")
    print(f"N2 (at {stress_ratio_2:.2f}*sigma_u): {N2:.3f} cycles")
    print(f"N3 (at {stress_ratio_3:.2f}*sigma_u): {N3:.3f} cycles\n")

    # Print the final equation used
    print("Palmgren-Miner Rule Calculation:")
    print(f"Fatigue Life = 1 / (({frac1} / {N1:.3f}) + ({frac2} / {N2:.3f}) + ({frac3} / {N3:.3f}))")
    print(f"Fatigue Life = 1 / ({damage_sum:.6f})")

    # Final result
    print("\n--- Final Answer ---")
    print(f"Estimated total fatigue life: {total_fatigue_life:.3f} cycles")
    
    # Format for the required Aeb scientific notation output
    final_answer_formatted = "{:.3e}".format(total_fatigue_life)
    print(f"In scientific notation (Aeb form): {final_answer_formatted}")
    print(f"\n<<<{final_answer_formatted}>>>")


if __name__ == "__main__":
    estimate_fatigue_life()