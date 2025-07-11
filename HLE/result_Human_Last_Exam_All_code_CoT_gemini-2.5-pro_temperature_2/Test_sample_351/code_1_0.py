def calculate_maintenance_fluid_needs():
    """
    Calculates maintenance fluid needs for a pediatric patient using the
    Body Surface Area (BSA) method and prints the steps.
    """
    # Patient and formula parameters
    bsa = 0.95  # Patient's Body Surface Area in m^2
    fluid_rate_per_bsa = 1500  # Standard rate in mL/m^2/day
    hours_in_day = 24

    # Step 1: Calculate total daily fluid volume
    daily_fluid_volume = fluid_rate_per_bsa * bsa

    # Step 2: Calculate the hourly fluid rate
    hourly_fluid_rate = daily_fluid_volume / hours_in_day

    print("The maintenance fluid rate is calculated using the BSA method (1500 mL/m^2/day).\n")

    # Print the equation for daily fluid calculation
    print("Step 1: Calculate the total fluid volume per day.")
    print(f"   {fluid_rate_per_bsa} mL/m^2 * {bsa} m^2 = {daily_fluid_volume:.2f} mL/day\n")

    # Print the equation for hourly fluid calculation
    print("Step 2: Convert the daily volume to an hourly rate.")
    print(f"   {daily_fluid_volume:.2f} mL / {hours_in_day} hours = {hourly_fluid_rate:.2f} cc/hr")


calculate_maintenance_fluid_needs()
<<<59.38>>>