def calculate_maintenance_fluids():
    """
    Calculates the maintenance fluid needs for a pediatric burn patient
    using the Body Surface Area (BSA) method.
    """
    # Patient data
    bsa = 0.95  # Body Surface Area in m^2

    # Standard calculation constants
    fluid_rate_per_m2_per_day = 1500  # mL/m^2/day
    hours_in_a_day = 24

    # 1. Calculate total daily fluid needs
    total_daily_fluids = fluid_rate_per_m2_per_day * bsa

    # 2. Calculate the hourly rate
    hourly_rate = total_daily_fluids / hours_in_a_day

    # Print the explanation and the final equation
    print("The maintenance fluid requirement is calculated using the formula: 1500 mL/m^2/day.")
    print(f"\nPatient's Body Surface Area (BSA): {bsa} m^2")
    print("\nCalculation:")
    print(f"({fluid_rate_per_m2_per_day} mL/m^2/day * {bsa} m^2) / {hours_in_a_day} hours")
    
    # Print the final result
    print(f"\nThe required maintenance fluid rate is {hourly_rate:.2f} cc/hr.")

calculate_maintenance_fluids()