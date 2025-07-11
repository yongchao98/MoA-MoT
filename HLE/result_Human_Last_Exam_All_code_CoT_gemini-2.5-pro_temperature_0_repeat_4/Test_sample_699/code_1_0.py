import math

def estimate_fatigue_life():
    """
    Estimates fatigue life using the Basquin relationship and Palmgren-Miner rule.
    """
    # Given parameters
    b = -0.09  # Basquin exponent

    # Life fractions for each stress level
    life_fractions = [0.70, 0.20, 0.10]

    # Stress levels as multiples of the endurance limit (sigma_e)
    stress_multiples = [1.0, 1.1, 1.2]

    # Assumption: sigma_e = 0.5 * sigma_uts (a common approximation for steels)
    # This allows us to express stress levels as a ratio of sigma_uts
    stress_ratios = [multiple * 0.5 for multiple in stress_multiples]

    # Calculate Nf for each stress level using the derived formula:
    # N_f = 0.5 * (sigma_a / sigma_uts)^(1/b)
    N_f_values = [0.5 * (ratio)**(1/b) for ratio in stress_ratios]
    
    N_f1, N_f2, N_f3 = N_f_values

    # Calculate the damage contribution from each stress level per total life L
    # Damage per L = (life_fraction_1 / N_f1) + (life_fraction_2 / N_f2) + ...
    damage_contributions = [lf / nf for lf, nf in zip(life_fractions, N_f_values)]
    
    d1, d2, d3 = damage_contributions

    # Total damage per L is the sum of individual contributions
    total_damage_per_L = sum(damage_contributions)

    # Total life L is the reciprocal of the total damage per L
    L = 1 / total_damage_per_L

    # --- Outputting the results step-by-step ---
    print("Step 1: Calculate cycles to failure (N_f) for each stress level.")
    print(f"   - At 1.0*σ_e (stress ratio = {stress_ratios[0]:.2f}): N_f1 = {N_f1:.2f} cycles")
    print(f"   - At 1.1*σ_e (stress ratio = {stress_ratios[1]:.2f}): N_f2 = {N_f2:.2f} cycles")
    print(f"   - At 1.2*σ_e (stress ratio = {stress_ratios[2]:.2f}): N_f3 = {N_f3:.2f} cycles")
    print("\nStep 2: Apply Palmgren-Miner rule: Σ (n_i / N_fi) = 1")
    print("   L * ( (life_fraction_1 / N_f1) + (life_fraction_2 / N_f2) + (life_fraction_3 / N_f3) ) = 1")
    print("   L = 1 / ( (0.70 / N_f1) + (0.20 / N_f2) + (0.10 / N_f3) )")
    print("\nStep 3: Substitute values and solve for L.")
    print(f"   L = 1 / ( (0.70 / {N_f1:.2f}) + (0.20 / {N_f2:.2f}) + (0.10 / {N_f3:.2f}) )")
    print(f"   L = 1 / ( {d1:.6f} + {d2:.6f} + {d3:.6f} )")
    print(f"   L = 1 / ( {total_damage_per_L:.6f} )")
    print(f"   L = {L:.2f} cycles")

    # Format the final answer in scientific notation Aeb (e.g., 2.665e+02)
    final_answer_formatted = "{:.3e}".format(L)
    print(f"\nEstimated Fatigue Life (L) in scientific notation: {final_answer_formatted} cycles")

# Execute the function
estimate_fatigue_life()
<<<2.665e+02>>>