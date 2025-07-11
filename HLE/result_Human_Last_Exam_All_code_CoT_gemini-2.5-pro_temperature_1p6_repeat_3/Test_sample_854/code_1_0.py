def calculate_antibiotic_degradation():
    """
    Calculates the remaining percentage of a heat-labile antibiotic after
    autoclaving, based on first-order decay kinetics.
    """
    # --- Input values from the problem description ---
    autoclave_time_minutes = 25
    
    # --- Assumptions for the model ---
    # The initial amount of antibiotic is represented as 100%.
    initial_concentration_percentage = 100.0
    # The half-life of chloramphenicol at 121°C is estimated to be 5 minutes.
    # This is a conservative estimate to demonstrate the degradation effect.
    half_life_at_121C_minutes = 5
    
    # --- Calculation Steps ---
    # The formula for half-life decay is: Final = Initial * (0.5)^(time / half-life)
    
    # 1. Calculate the number of half-lives that occur during the autoclave cycle.
    num_half_lives = autoclave_time_minutes / half_life_at_121C_minutes
    
    # 2. Calculate the final concentration percentage.
    final_concentration_percentage = initial_concentration_percentage * (0.5) ** num_half_lives
    
    # --- Output the results ---
    print("--- Analysis of Antibiotic Inactivation in Batch 3 ---")
    print(f"The PDA was autoclaved for {autoclave_time_minutes} minutes.")
    print(f"The estimated half-life of the antibiotic at autoclave temperature is {half_life_at_121C_minutes} minutes.")
    print("\nThe calculation for the remaining antibiotic is based on the decay formula: Final = Initial * (0.5)^(time / half_life)")
    print(f"Final % = {initial_concentration_percentage} * (0.5)^({autoclave_time_minutes} / {half_life_at_121C_minutes})")
    print(f"Final % = {initial_concentration_percentage} * (0.5)^({num_half_lives})")
    print(f"Final % = {final_concentration_percentage:.4f}%")
    
    print("\nConclusion:")
    print(f"After {autoclave_time_minutes} minutes of autoclaving, only {final_concentration_percentage:.4f}% of the antibiotic remains active.")
    print("This negligible amount is insufficient to inhibit bacterial growth, explaining the contamination in Batch 3.")

calculate_antibiotic_degradation()