import math

def solve_microbiology_case():
    """
    Analyzes the microbiology lab scenario by simulating bacterial growth.
    """
    # Assign arbitrary growth rates. Contaminants often grow faster than pathogens like Campylobacter.
    # Growth rate is in arbitrary "size units" per hour.
    contaminant_growth_rate = 0.1  # e.g., Bacillus species
    pathogen_growth_rate = 0.06  # e.g., Campylobacter species

    # The colony must reach a certain size to be easily visible.
    visibility_threshold = 3.0

    # Initial observation time by the first lab.
    initial_incubation_hours = 48

    # --- Analysis at the initial 48-hour observation ---
    print("--- Analysis at 48 Hours ---")
    contaminant_size_48h = contaminant_growth_rate * initial_incubation_hours
    pathogen_size_48h = pathogen_growth_rate * initial_incubation_hours

    print(f"Equation for Contaminant Size at {initial_incubation_hours} hours: {contaminant_growth_rate} * {initial_incubation_hours} = {contaminant_size_48h:.2f}")
    print(f"Equation for Pathogen Size at {initial_incubation_hours} hours: {pathogen_growth_rate} * {initial_incubation_hours} = {pathogen_size_48h:.2f}")

    if contaminant_size_48h > visibility_threshold:
        print(f"Result: The contaminant colony is visible (size {contaminant_size_48h:.2f} > threshold {visibility_threshold}).")
    else:
        print(f"Result: The contaminant colony is not yet visible.")

    if pathogen_size_48h > visibility_threshold:
        print(f"Result: The pathogen colony is visible.")
    else:
        print(f"Result: The pathogen colony is NOT yet visible (size {pathogen_size_48h:.2f} < threshold {visibility_threshold}).")

    print("\nConclusion at 48 hours: The lab likely saw the fast-growing contaminant but missed the slower-growing pathogen.")
    
    # --- Analysis with extended incubation ---
    # Let's see what happens if they incubate for one more day.
    longer_incubation_hours = 72

    print(f"\n--- Analysis at {longer_incubation_hours} Hours (Extended Incubation) ---")
    contaminant_size_72h = contaminant_growth_rate * longer_incubation_hours
    pathogen_size_72h = pathogen_growth_rate * longer_incubation_hours
    
    print(f"Equation for Contaminant Size at {longer_incubation_hours} hours: {contaminant_growth_rate} * {longer_incubation_hours} = {contaminant_size_72h:.2f}")
    print(f"Equation for Pathogen Size at {longer_incubation_hours} hours: {pathogen_growth_rate} * {longer_incubation_hours} = {pathogen_size_72h:.2f}")

    if pathogen_size_72h > visibility_threshold:
        print(f"Result: After extended incubation, the pathogen colony is now visible (size {pathogen_size_72h:.2f} > threshold {visibility_threshold}).")
    else:
        print("Result: Even after extended incubation, the pathogen is not visible.")

    print("\nFinal Conclusion: By incubating the sample for longer, the slower-growing Campylobacter would have had time to form visible colonies.")


solve_microbiology_case()