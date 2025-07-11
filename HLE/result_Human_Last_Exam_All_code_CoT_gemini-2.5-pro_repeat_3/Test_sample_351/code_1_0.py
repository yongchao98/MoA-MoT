def calculate_maintenance_fluids():
    """
    Calculates pediatric maintenance fluid needs using the Body Surface Area (BSA) method.
    """
    # Patient data
    bsa = 0.95  # Body Surface Area in m^2
    
    # Constants for calculation
    fluid_rate_per_m2 = 1500  # Standard rate in mL/m^2/day
    hours_in_day = 24         # Hours in a day

    # Calculate total daily fluid volume
    total_daily_fluid = bsa * fluid_rate_per_m2

    # Calculate the hourly rate
    hourly_rate = total_daily_fluid / hours_in_day

    # Print the explanation and the final calculation
    print("Calculating maintenance fluid needs using the Body Surface Area (BSA) method.")
    print("The formula is: (BSA * 1500 mL/m^2) / 24 hours\n")
    print("Final Calculation:")
    # The final equation with each number explicitly shown
    print(f"({bsa} m^2 * {fluid_rate_per_m2} mL/m^2) / {hours_in_day} hours = {hourly_rate:.3f} cc/hr")

# Execute the function
calculate_maintenance_fluids()