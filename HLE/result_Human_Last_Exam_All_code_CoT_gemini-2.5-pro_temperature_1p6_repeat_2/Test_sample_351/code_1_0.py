def calculate_maintenance_fluids():
    """
    Calculates maintenance fluid needs for a pediatric patient using the BSA method.
    """
    # Patient data
    bsa = 0.95  # Body Surface Area in m^2
    
    # Constants
    fluid_per_m2_per_day = 1500  # Standard constant in ml/m^2/day
    hours_in_day = 24
    
    # Step 1: Calculate total daily maintenance fluid
    daily_fluid = fluid_per_m2_per_day * bsa
    
    # Step 2: Calculate hourly maintenance fluid
    hourly_fluid = daily_fluid / hours_in_day
    
    # Output the calculation steps and the final answer
    print(f"Calculating daily maintenance fluid based on BSA:")
    print(f"Equation: Fluid (ml/day) = Constant (ml/m^2/day) * BSA (m^2)")
    print(f"Calculation: Daily Fluid = {fluid_per_m2_per_day} * {bsa} = {daily_fluid:.2f} ml/day")
    print("\nCalculating hourly maintenance fluid rate:")
    print(f"Equation: Fluid (cc/hr) = Daily Fluid (ml) / 24 hours")
    print(f"Calculation: Hourly Fluid = {daily_fluid:.2f} / {hours_in_day} = {hourly_fluid:.2f} cc/hr")

calculate_maintenance_fluids()