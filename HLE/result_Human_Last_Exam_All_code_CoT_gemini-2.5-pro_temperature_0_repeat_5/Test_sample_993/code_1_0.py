def analyze_paleoclimate_conditions():
    """
    Analyzes environmental factors to determine which would lead to low levels
    of the 18O isotope in lake sulfate deposits.
    """

    print("Goal: Identify conditions leading to low 18O levels in ancient lake sulfate.")
    print("="*70)
    print("Key Scientific Principles:")
    print("1. Evaporation: Lighter 16O evaporates more easily than heavier 18O. High evaporation enriches the remaining lake water in 18O.")
    print("2. Inflow: Precipitation (rain/snow) is naturally low in 18O. High inflow dilutes the lake water, keeping 18O levels low.")
    print("3. Sulfate Formation: Sulfate minerals record the 18O/16O ratio of the lake water they form in.")
    print("\nConclusion: To get low 18O in sulfate, we need low 18O in the lake water. This requires minimizing evaporation and maximizing inflow.")
    print("="*70)

    # Define how each factor influences 18O levels
    # A lower score means more favorable for low 18O
    factor_effects = {
        'Climate': {'Wet': 'Favorable (reduces evaporation, increases inflow)', 'Dry': 'Unfavorable (increases evaporation)'},
        'Temperature': {'Cold': 'Favorable (reduces evaporation)', 'Warm': 'Unfavorable (increases evaporation)'},
        'Lake Level': {'High': 'Favorable (buffers against evaporation effects)', 'Shallow': 'Unfavorable (sensitive to evaporation)'}
    }

    # All possible choices from the study
    options = {
        'A': ('Wet', 'warm', 'shallow'),
        'B': ('Dry', 'warm', 'shallow'),
        'C': ('Wet', 'cold', 'shallow'),
        'D': ('Dry', 'cold', 'shallow'),
        'E': ('Wet', 'warm', 'high'),
        'F': ('Dry', 'warm', 'high'),
        'G': ('Wet', 'cold', 'high'),
        'H': ('Dry', 'cold', 'high')
    }

    print("Evaluating the choices:\n")
    best_choice = None
    best_score = float('inf')

    for key, (climate, temp, level) in options.items():
        score = 0
        if climate.lower() == 'dry': score += 1
        if temp.lower() == 'warm': score += 1
        if level.lower() == 'shallow': score += 1
        
        if score < best_score:
            best_score = score
            best_choice = key

    print(f"The conditions most favorable for low 18O are: Wet, Cold, and High lake levels.")
    print(f"This corresponds to the option with the most 'Favorable' ratings.\n")
    print(f"The best choice is '{best_choice}', which represents: {options[best_choice][0]}, {options[best_choice][1]} climate with {options[best_choice][2]} lake levels.")

analyze_paleoclimate_conditions()