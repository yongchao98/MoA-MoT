def solve():
    """
    Models and explains the solution to the laboratory scenario.
    """
    print("Analyzing the laboratory scenario by modeling bacterial growth...")
    print("="*60)

    # --- Assumptions for the Model ---
    # We will model exponential growth based on a doubling time.
    # These values are for demonstration purposes.
    contaminant_doubling_time_h = 8  # A fast-growing organism
    pathogen_doubling_time_h = 12 # Campylobacter is a slower-growing organism

    # --- Observation Time Points ---
    initial_incubation_h = 48
    extended_incubation_h = 72

    # --- The Growth Equation ---
    # final_size = initial_size * (2^(time / doubling_time))
    # We'll assume an initial_size of 1 for simplicity.
    print(f"The model uses the equation: Final Size = 1 * (2 ^ (Time / Doubling Time))")
    print(f"Contaminant (e.g., Bacillus) doubles every {contaminant_doubling_time_h} hours.")
    print(f"Pathogen (Campylobacter) doubles every {pathogen_doubling_time_h} hours.")
    print("-" * 60)

    # --- Calculation at 48 Hours ---
    print(f"CALCULATING COLONY SIZE AT {initial_incubation_h} HOURS:")
    # Contaminant
    contaminant_size_48h = 2**(initial_incubation_h / contaminant_doubling_time_h)
    print(f"  Contaminant Size = 2 ^ ({initial_incubation_h} / {contaminant_doubling_time_h}) = {contaminant_size_48h:.1f}")

    # Pathogen
    pathogen_size_48h = 2**(initial_incubation_h / pathogen_doubling_time_h)
    print(f"  Pathogen Size    = 2 ^ ({initial_incubation_h} / {pathogen_doubling_time_h}) = {pathogen_size_48h:.1f}")
    print("\nObservation: At 48 hours, the contaminant colony is significantly larger than the pathogen's. The small pathogen colonies could be easily missed or obscured.")
    print("-" * 60)


    # --- Calculation at 72 Hours ---
    print(f"CALCULATING COLONY SIZE AT {extended_incubation_h} HOURS:")
    # Contaminant
    contaminant_size_72h = 2**(extended_incubation_h / contaminant_doubling_time_h)
    print(f"  Contaminant Size = 2 ^ ({extended_incubation_h} / {contaminant_doubling_time_h}) = {contaminant_size_72h:.1f}")

    # Pathogen
    pathogen_size_72h = 2**(extended_incubation_h / pathogen_doubling_time_h)
    print(f"  Pathogen Size    = 2 ^ ({extended_incubation_h} / {pathogen_doubling_time_h}) = {pathogen_size_72h:.1f}")
    print("\nObservation: After extending incubation, the pathogen's colony has grown to a much more noticeable size, increasing the likelihood of being seen and identified correctly.")
    print("="*60)
    print("Conclusion: The simulation shows that incubating the plates for a longer duration is the most effective way to recover a slow-growing pathogen that was initially missed.")


solve()
<<<D>>>