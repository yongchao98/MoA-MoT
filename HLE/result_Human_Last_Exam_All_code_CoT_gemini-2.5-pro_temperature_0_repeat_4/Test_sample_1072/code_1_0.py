def calculate_ages():
    """
    Calculates the (U-Th)/He ages for the three geological samples based on their thermal histories.
    """
    T_surface = 25  # °C
    geotherm = 25  # °C/km
    Tc_zircon = 180  # °C, closure temperature for Zircon (U-Th)/He
    Tc_apatite = 75  # °C, closure temperature for Apatite (U-Th)/He

    print("--- Analysis of Sample Ages ---")

    # --- Sample 1: Zircon from a pluton ---
    print("\n[Sample 1: Zircon]")
    depth_1_start = 15  # km
    time_1_start = 100  # Ma
    
    # Temperature at the start of exhumation
    T_1_start = (depth_1_start * geotherm) + T_surface
    print(f"Initial temperature at {depth_1_start} km depth at {time_1_start} Ma = ({depth_1_start} km * {geotherm}°C/km) + {T_surface}°C = {T_1_start}°C")
    
    # Exhumation rate
    exhumation_rate_1 = depth_1_start / time_1_start  # km/Myr
    print(f"Exhumation rate = {depth_1_start} km / {time_1_start} Myr = {exhumation_rate_1} km/Myr")

    # Cooling rate
    cooling_rate_1 = exhumation_rate_1 * geotherm # °C/Myr
    print(f"Cooling rate = {exhumation_rate_1:.2f} km/Myr * {geotherm}°C/km = {cooling_rate_1}°C/Myr")

    # Time to cool from initial temperature to closure temperature
    time_to_cool_1 = (T_1_start - Tc_zircon) / cooling_rate_1
    print(f"Time to cool from {T_1_start}°C to {Tc_zircon}°C = ({T_1_start} - {Tc_zircon})°C / {cooling_rate_1}°C/Myr = {time_to_cool_1:.1f} Myr")

    # Final (U-Th)/He date
    age_1 = time_1_start - time_to_cool_1
    print(f"Calculated Age = {time_1_start} Ma - {time_to_cool_1:.1f} Myr = {age_1:.1f} Ma")

    # --- Sample 2: Apatite from a sedimentary rock ---
    print("\n[Sample 2: Apatite]")
    T_2_start = 250 # °C
    time_2_start = 100 # Ma
    
    # The clock is reset at 100 Ma at 250°C. First, find the depth.
    depth_2_start = (T_2_start - T_surface) / geotherm
    print(f"Initial depth at {time_2_start} Ma = ({T_2_start}°C - {T_surface}°C) / {geotherm}°C/km = {depth_2_start} km")

    # Exhumation rate
    exhumation_rate_2 = depth_2_start / time_2_start # km/Myr
    print(f"Exhumation rate = {depth_2_start} km / {time_2_start} Myr = {exhumation_rate_2} km/Myr")

    # Cooling rate
    cooling_rate_2 = exhumation_rate_2 * geotherm # °C/Myr
    print(f"Cooling rate = {exhumation_rate_2:.2f} km/Myr * {geotherm}°C/km = {cooling_rate_2}°C/Myr")

    # Time to cool from initial temperature to closure temperature
    time_to_cool_2 = (T_2_start - Tc_apatite) / cooling_rate_2
    print(f"Time to cool from {T_2_start}°C to {Tc_apatite}°C = ({T_2_start} - {Tc_apatite})°C / {cooling_rate_2}°C/Myr = {time_to_cool_2:.1f} Myr")

    # Final (U-Th)/He date
    age_2 = time_2_start - time_to_cool_2
    print(f"Calculated Age = {time_2_start} Ma - {time_to_cool_2:.1f} Myr = {age_2:.1f} Ma")

    # --- Sample 3: Apatite from a rhyolite ---
    print("\n[Sample 3: Apatite]")
    age_3 = 90  # Ma (eruption age)
    print(f"Age from eruption at 90 Ma = {age_3} Ma")

    # --- Conclusion on Ages ---
    print("\n--- Age Comparison ---")
    print(f"Sample 1 Age: ~{age_1:.1f} Ma")
    print(f"Sample 2 Age: ~{age_2:.1f} Ma")
    print(f"Sample 3 Age: {age_3} Ma")
    print("\nConclusion: Sample 3 is the oldest, and Sample 2 is the youngest.")
    print("Therefore, statement [H] is TRUE.")

calculate_ages()