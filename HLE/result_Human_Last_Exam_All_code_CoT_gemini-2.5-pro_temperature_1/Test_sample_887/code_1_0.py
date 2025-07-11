import math

def analyze_reaction_selectivity():
    """
    This script models the effect of vibrational excitation on the reaction
    selectivity in F + CHD3, demonstrating bond-selective chemistry.
    """

    # --- Model Parameters (Illustrative Values) ---
    # Baseline activation energy for C-H cleavage (kJ/mol)
    Ea_CH_base = 5.0
    # Baseline activation energy for C-D cleavage (kJ/mol). C-D bond is stronger.
    Ea_CD_base = 6.0
    # Vibrational energy added to the C-H bond by the laser (kJ/mol)
    E_vib_laser = 4.5
    # Thermal energy factor (RT in kJ/mol, corresponds to ~300 K)
    RT = 2.5

    print("This model demonstrates how exciting the C-H bond affects its reactivity.")
    print("-" * 60)

    # --- Case 1: No Laser Excitation (Thermal Reaction) ---
    print("Scenario 1: No Laser Excitation\n")
    # Rate factors are proportional to exp(-Ea/RT) from the Arrhenius equation.
    k_factor_CH_no_laser = math.exp(-Ea_CH_base / RT)
    k_factor_CD_no_laser = math.exp(-Ea_CD_base / RT)
    
    # The ratio of rates determines the product branching ratio (H removal vs. D removal)
    if k_factor_CD_no_laser > 0:
        ratio_no_laser = k_factor_CH_no_laser / k_factor_CD_no_laser
        print(f"Activation Energy (C-H): {Ea_CH_base:.1f} kJ/mol")
        print(f"Activation Energy (C-D): {Ea_CD_base:.1f} kJ/mol")
        print(f"Initial relative rate of H vs. D removal: {ratio_no_laser:.2f}")
        print("Without the laser, H removal is slightly favored due to its lower activation energy.")
    else:
        print("C-D reaction rate is negligible.")
    
    print("-" * 60)

    # --- Case 2: With Laser Excitation of C-H bond ---
    print("Scenario 2: With Laser Excitation of the C-H Bond\n")
    # The laser energy is assumed to directly lower the C-H activation barrier.
    Ea_CH_effective = max(0, Ea_CH_base - E_vib_laser)
    # The C-D bond activation energy is unaffected.
    Ea_CD_effective = Ea_CD_base

    k_factor_CH_with_laser = math.exp(-Ea_CH_effective / RT)
    k_factor_CD_with_laser = math.exp(-Ea_CD_effective / RT) # This is unchanged

    # Final Equation: New Ratio = exp(-(Ea_CH-E_vib)/RT) / exp(-Ea_CD/RT)
    print("Calculating the new relative rate:")
    print(f"New relative rate = exp(-({Ea_CH_base} - {E_vib_laser}) / {RT}) / exp(-{Ea_CD_effective} / {RT})")

    if k_factor_CD_with_laser > 0:
        ratio_with_laser = k_factor_CH_with_laser / k_factor_CD_with_laser
        print(f"\nEffective Activation Energy (C-H): {Ea_CH_effective:.1f} kJ/mol")
        print(f"Effective Activation Energy (C-D): {Ea_CD_effective:.1f} kJ/mol (unchanged)")
        print(f"New relative rate of H vs. D removal: {ratio_with_laser:.2f}")
    
        enhancement_factor = ratio_with_laser / ratio_no_laser
        print(f"\nThe laser enhances the preference for H atom removal by a factor of {enhancement_factor:.1f}.")
    else:
        print("C-D reaction rate is negligible.")

    print("\nConclusion: Exciting the C-H bond accelerates the reaction by specifically")
    print("enhancing the likelihood of H atom removal over D atoms.")

analyze_reaction_selectivity()