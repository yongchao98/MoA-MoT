import math

def estimate_fatigue_life():
    """
    Estimates fatigue life using the Palmgren-Miner rule based on the given
    Basquin relationship and loading conditions.
    """
    # --- Given parameters and assumptions ---
    b = -0.09  # Basquin exponent

    # Cycle fractions for the three stress levels
    c1_fraction = 0.70
    c2_fraction = 0.20
    c3_fraction = 0.10

    # Stress levels as a ratio of the endurance limit (σ_e)
    sigma_ratio_2 = 1.1
    sigma_ratio_3 = 1.2

    # Assumption: The endurance limit (σ_e) corresponds to a life of 1,000,000 cycles.
    # Stresses at or below σ_e are considered non-damaging.
    N_endurance = 1e6

    # --- Calculations ---

    # The life at the endurance limit (N_f1) is effectively infinite for damage calculation.
    # The damage contribution from the first stress level is zero.

    # Calculate the number of cycles to failure (N_f) for the higher stress levels
    # using the relative form of the Basquin equation: N_f = N_e * (σ/σ_e)^(1/b)
    inv_b = 1.0 / b
    N_f2 = N_endurance * (sigma_ratio_2 ** inv_b)
    N_f3 = N_endurance * (sigma_ratio_3 ** inv_b)

    # Apply the Palmgren-Miner rule to find the total life (N_total).
    # The damage equation is: 1 = N_total * (c1/N_f1 + c2/N_f2 + c3/N_f3)
    # Since the first term is zero: 1 = N_total * (c2/N_f2 + c3/N_f3)
    damage_sum_per_cycle = (c2_fraction / N_f2) + (c3_fraction / N_f3)
    N_total_life = 1 / damage_sum_per_cycle

    # --- Output ---
    print("--- Fatigue Life Estimation ---")
    print(f"Assumed life at endurance limit (N_e): {N_endurance:.0e} cycles")
    print(f"Life at 1.1*σ_e (N_f2): {N_f2:.3e} cycles")
    print(f"Life at 1.2*σ_e (N_f3): {N_f3:.3e} cycles")
    print("-" * 35)

    # Print the final equation with all numbers, as requested.
    print("Final Palmgren-Miner Damage Equation:")
    # The equation is D = 1 = N_total * (c1/N_f1 + c2/N_f2 + c3/N_f3)
    print(f"1 = N_total * ( ({c1_fraction} / ∞) + ({c2_fraction} / {N_f2:.3e}) + ({c3_fraction} / {N_f3:.3e}) )")
    print("-" * 35)

    # Format the final answer to 3 decimal places in scientific notation (Aeb form)
    formatted_life = "{:.3e}".format(N_total_life)
    print(f"Estimated Total Fatigue Life (N_total): {formatted_life} cycles")

    # Return the final answer in the specified format
    print(f"\n<<<{formatted_life}>>>")

# Execute the function
estimate_fatigue_life()