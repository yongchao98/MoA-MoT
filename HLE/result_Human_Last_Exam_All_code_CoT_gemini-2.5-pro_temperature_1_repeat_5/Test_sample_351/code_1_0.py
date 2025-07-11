def calculate_maintenance_fluids():
    """
    Calculates maintenance fluid needs for a pediatric patient based on Body Surface Area (BSA).
    """
    # Patient Information
    bsa = 0.95  # Body Surface Area in m^2
    
    # Constants
    fluid_rate_per_m2_per_day = 1500  # Standard maintenance fluid rate in mL/m^2/day
    hours_in_day = 24

    # Step 1: Calculate total daily fluid needs
    daily_fluid_needs = fluid_rate_per_m2_per_day * bsa

    # Step 2: Calculate hourly fluid needs
    hourly_fluid_needs = daily_fluid_needs / hours_in_day

    # Step 3: Print the results
    print("To calculate maintenance fluid needs using the Body Surface Area (BSA) method:")
    print("Formula: (1500 mL/m^2/day * BSA) / 24 hours/day\n")
    
    print("1. Calculate Daily Fluid Needs:")
    print(f"{fluid_rate_per_m2_per_day} mL/m^2/day * {bsa} m^2 = {daily_fluid_needs:.2f} mL/day\n")

    print("2. Calculate Hourly Fluid Needs:")
    print(f"{daily_fluid_needs:.2f} mL/day / {hours_in_day} hours = {hourly_fluid_needs:.2f} cc/hr\n")

    print(f"The calculated maintenance fluid rate is {hourly_fluid_needs:.2f} cc/hr.")

if __name__ == "__main__":
    calculate_maintenance_fluids()