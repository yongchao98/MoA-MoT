import math

def estimate_fatigue_life():
    """
    Estimates fatigue life using Basquin's relation and Palmgren-Miner rule.
    """
    # Given constants and loading fractions
    b = -0.09  # Basquin exponent
    cycle_fractions = [0.70, 0.20, 0.10]
    stress_multipliers = [1.0, 1.1, 1.2] # Stress levels as multipliers of sigma_e

    # --- Step 1 & 2: Derive S-N parameters and make an assumption ---
    # The S-N relationship is σ_a = C * N^b.
    # We are given that σ_a = σ_uts when N = 0.5.
    # This leads to a general formula for life N at a given stress σ_a:
    # N(σ_a) = 0.5 * (σ_a / σ_uts)^(1/b)
    #
    # To solve, we must assume a relationship between the endurance limit (σ_e)
    # and the ultimate tensile strength (σ_uts). A common assumption for steels is used.
    sigma_e_to_uts_ratio = 0.5
    print(f"Plan:")
    print(f"1. Use Basquin's relationship: σ_a = C * N^b, with b = {b}")
    print(f"2. Assume endurance limit σ_e = {sigma_e_to_uts_ratio} * σ_uts")
    print(f"3. Calculate life N_i for each stress level σ_i using N_i = 0.5 * (σ_i / σ_uts)^(1/b)")
    print(f"4. Use Palmgren-Miner rule to find total life: N_total = 1 / Σ(n_i / N_i)\n")


    # --- Step 3: Calculate life (N_f) for each stress level ---
    # The stress ratios (σ_i / σ_uts) are calculated as: multiplier * (σ_e / σ_uts)
    stress_ratios = [m * sigma_e_to_uts_ratio for m in stress_multipliers]

    def calculate_Nf(stress_ratio, b_exp):
        """Calculates cycles to failure (Nf) based on a stress ratio relative to UTS."""
        return 0.5 * (stress_ratio ** (1 / b_exp))

    N1 = calculate_Nf(stress_ratios[0], b) # Life at 1.0*σ_e
    N2 = calculate_Nf(stress_ratios[1], b) # Life at 1.1*σ_e
    N3 = calculate_Nf(stress_ratios[2], b) # Life at 1.2*σ_e

    print("--- Calculation ---")
    print(f"Cycles to failure for each stress level:")
    print(f"N1 (at {stress_multipliers[0]:.1f}*σ_e): {N1:.3f} cycles")
    print(f"N2 (at {stress_multipliers[1]:.1f}*σ_e): {N2:.3f} cycles")
    print(f"N3 (at {stress_multipliers[2]:.1f}*σ_e): {N3:.3f} cycles\n")

    # --- Step 4: Apply Palmgren-Miner Rule ---
    # The total damage D = 1 = (n1/N1) + (n2/N2) + (n3/N3)
    # Since n_i are fractions of total life N_total (e.g., n1 = 0.70 * N_total),
    # we get N_total = 1 / (0.70/N1 + 0.20/N2 + 0.10/N3)
    damage_sum_per_block = (cycle_fractions[0] / N1) + \
                           (cycle_fractions[1] / N2) + \
                           (cycle_fractions[2] / N3)

    total_life_Nf = 1 / damage_sum_per_block

    # --- Step 5: Format and Print Output ---
    print("Final Equation using Palmgren-Miner Rule:")
    equation_str = (
        f"N_total = 1 / ( {cycle_fractions[0]:.2f} / {N1:.3f} + "
        f"{cycle_fractions[1]:.2f} / {N2:.3f} + "
        f"{cycle_fractions[2]:.2f} / {N3:.3f} )"
    )
    print(equation_str)
    
    result_str = f"N_total = {total_life_Nf:.3f} cycles"
    print(result_str)

    # Final answer in scientific notation
    final_answer_sci = f"{total_life_Nf:.3e}"
    print(f"\nEstimated Total Fatigue Life (Scientific Notation): {final_answer_sci}")


estimate_fatigue_life()
<<<6.289e+02>>>