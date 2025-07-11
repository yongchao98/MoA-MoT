import math

def estimate_fatigue_life():
    """
    Estimates fatigue life using Basquin's law and Palmgren-Miner rule.
    """
    # Given parameters
    b = -0.09  # Basquin exponent

    # Cycle fractions for each stress level
    frac_1 = 0.70  # 70% of life
    frac_2 = 0.20  # 20% of life
    frac_3 = 0.10  # 10% of life

    # --- Step 1: Define stress ratios based on the assumption that sigma_e = 0.5 * sigma_uts ---
    # We don't need the absolute values of sigma_e or sigma_uts, only their ratio.
    # Stress Ratio = sigma_a / sigma_uts
    # Assumption: sigma_e (endurance limit) = 0.5 * sigma_uts
    
    # Stress level 1: sigma_e
    stress_ratio_1 = 0.50
    # Stress level 2: 1.1 * sigma_e
    stress_ratio_2 = 1.1 * 0.50
    # Stress level 3: 1.2 * sigma_e
    stress_ratio_3 = 1.2 * 0.50

    # --- Step 2: Calculate cycles to failure (N) for each stress level ---
    # From Basquin's law and the condition that N=0.5 when sigma_a = sigma_uts,
    # we derive N_f = 0.5 * (sigma_a / sigma_uts)^(1/b)
    
    N1 = 0.5 * (stress_ratio_1)**(1 / b)
    N2 = 0.5 * (stress_ratio_2)**(1 / b)
    N3 = 0.5 * (stress_ratio_3)**(1 / b)

    # --- Step 3: Apply Palmgren-Miner linear damage rule ---
    # Total Life N_total = 1 / sum(fraction_i / N_i)
    total_life = 1 / ((frac_1 / N1) + (frac_2 / N2) + (frac_3 / N3))

    # --- Step 4: Output the results as requested ---
    # Print the equation with all the numbers
    print("This solution assumes the fatigue endurance limit is 50% of the ultimate tensile strength.")
    print("\nCalculation using the Palmgren-Miner rule:")
    print(f"N_total = 1 / (({frac_1}/{N1:.1f}) + ({frac_2}/{N2:.1f}) + ({frac_3}/{N3:.1f}))")
    
    # Print the final calculated life with formatting
    print(f"\nEstimated total fatigue life: {total_life:.3f} cycles")
    
    # Format the final answer in scientific notation with 3 decimal places
    final_answer_sci = f"{total_life:.3e}"
    
    print(f"In scientific notation: {final_answer_sci}")
    
    # Return the answer in the requested format
    print(f"\n<<<{final_answer_sci}>>>")

# Execute the function
estimate_fatigue_life()
