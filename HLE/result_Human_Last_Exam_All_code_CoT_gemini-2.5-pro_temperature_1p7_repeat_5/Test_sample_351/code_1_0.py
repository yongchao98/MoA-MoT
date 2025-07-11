def calculate_maintenance_fluids():
    """
    Calculates and prints the maintenance fluid needs for a patient based on BSA.
    """
    # Patient and formula specific data
    bsa_m2 = 0.95  # Body Surface Area in square meters
    fluid_rate_per_m2_per_day = 1500  # Standard rate in mL/m^2/day
    hours_per_day = 24

    # 1. Calculate the total daily fluid needs
    daily_fluid_needs_ml = fluid_rate_per_m2_per_day * bsa_m2

    # 2. Convert the daily need to an hourly rate
    hourly_fluid_needs_cc_hr = daily_fluid_needs_ml / hours_per_day

    # 3. Print the full equation and the result
    # The final answer is requested in cc/hr, which is equivalent to mL/hr.
    print("Calculation of maintenance fluid needs using the BSA method:")
    print(f"({fluid_rate_per_m2_per_day} mL/m^2/day * {bsa_m2} m^2) / {hours_per_day} hours/day = {round(hourly_fluid_needs_cc_hr, 2)} cc/hr")

# Execute the calculation
calculate_maintenance_fluids()