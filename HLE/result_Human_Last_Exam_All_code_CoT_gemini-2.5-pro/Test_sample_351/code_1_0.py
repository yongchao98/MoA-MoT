def calculate_maintenance_fluids():
    """
    Calculates pediatric maintenance fluid needs based on Body Surface Area (BSA).
    """
    # Patient and formula constants
    bsa = 0.95  # Body Surface Area in m^2
    daily_rate_per_m2 = 1500  # Standard maintenance rate in cc/m^2/day
    hours_in_day = 24

    # Step 1: Calculate total daily maintenance fluid
    total_daily_fluid = daily_rate_per_m2 * bsa

    # Step 2: Calculate the hourly rate
    hourly_rate = total_daily_fluid / hours_in_day

    # Output the explanation and the final equation
    print("Calculating pediatric maintenance fluid needs using the Body Surface Area (BSA) method.")
    print("The formula is: (Standard Rate * BSA) / Hours in a day\n")
    
    # As requested, printing the final equation with each number
    print(f"The final equation is: ({daily_rate_per_m2} cc/m^2 * {bsa} m^2) / {hours_in_day} hours")
    
    print(f"\nThe calculated hourly maintenance fluid rate is: {hourly_rate:.2f} cc/hr")

calculate_maintenance_fluids()