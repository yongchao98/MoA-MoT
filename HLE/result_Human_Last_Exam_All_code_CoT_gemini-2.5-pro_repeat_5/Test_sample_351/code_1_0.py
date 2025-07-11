def calculate_maintenance_fluids():
    """
    Calculates maintenance fluid needs for a pediatric patient using the Body Surface Area (BSA) method.
    """
    # Patient data
    bsa = 0.95  # Body Surface Area in m^2

    # Constants for the BSA method
    fluid_per_m2_per_day = 1500  # Standard rate in cc/m^2/day
    hours_in_day = 24

    # Step 1: Calculate the total daily maintenance fluid requirement.
    daily_fluid_volume = fluid_per_m2_per_day * bsa

    # Step 2: Convert the daily requirement to an hourly rate.
    hourly_fluid_rate = daily_fluid_volume / hours_in_day

    print("Calculating Maintenance Fluid Needs using the BSA Method (1500 cc/m^2/day)")
    print("-" * 60)
    
    print("Step 1: Calculate total daily fluid volume.")
    print(f"Equation: Daily Volume = Rate * BSA")
    print(f"Calculation: {fluid_per_m2_per_day} cc/m^2/day * {bsa} m^2 = {daily_fluid_volume:.2f} cc/day")
    
    print("\nStep 2: Convert daily volume to an hourly rate.")
    print(f"Equation: Hourly Rate = Daily Volume / 24 hours")
    print(f"Calculation: {daily_fluid_volume:.2f} cc/day / {hours_in_day} hours = {hourly_fluid_rate:.2f} cc/hr")
    
    print("\nFinal Answer:")
    print(f"The patient's maintenance fluid rate is {hourly_fluid_rate:.2f} cc/hr.")

calculate_maintenance_fluids()