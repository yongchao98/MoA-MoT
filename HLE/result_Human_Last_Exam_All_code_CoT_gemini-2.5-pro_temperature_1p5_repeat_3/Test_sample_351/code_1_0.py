def calculate_maintenance_fluids():
    """
    Calculates maintenance fluid needs for a pediatric patient using the
    Body Surface Area (BSA) method.
    """
    # Patient Data
    bsa_m2 = 0.95  # Patient's Body Surface Area in square meters

    # Standard Formula Constants
    fluid_per_bsa_per_day = 1500  # Standard rate in mL/m^2/day
    hours_in_day = 24  # Number of hours in a day

    # Step 1: Calculate total daily fluid volume
    total_daily_volume_ml = fluid_per_bsa_per_day * bsa_m2

    # Step 2: Calculate the hourly rate in cc/hr (1 mL = 1 cc)
    hourly_rate_cc_hr = total_daily_volume_ml / hours_in_day

    # Output the explanation, the full equation, and the result
    print("Calculating maintenance fluid rate using the BSA method (1500 mL/m^2/day):")
    print(f"({fluid_per_bsa_per_day} mL/m^2/day * {bsa_m2} m^2) / {hours_in_day} hours = {hourly_rate_cc_hr:.2f} cc/hr")

# Execute the function
calculate_maintenance_fluids()