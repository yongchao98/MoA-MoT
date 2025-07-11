def calculate_maintenance_fluids():
    """
    Calculates pediatric maintenance fluid needs based on Body Surface Area (BSA).
    """
    # Patient data
    bsa = 0.95  # Body Surface Area in m^2
    
    # Standard constants
    daily_fluid_rate_per_m2 = 1500  # Standard rate in mL/m^2/day
    hours_in_day = 24
    
    # 1. Calculate total daily fluid volume
    total_daily_fluid = daily_fluid_rate_per_m2 * bsa
    
    # 2. Calculate the hourly rate
    # 1 mL = 1 cc, so mL/hr is equivalent to cc/hr
    hourly_rate = total_daily_fluid / hours_in_day
    
    # 3. Print the equation and the final answer
    print("Calculation for maintenance fluid needs in cc/hr:")
    print(f"({daily_fluid_rate_per_m2} mL/m^2/day * {bsa} m^2) / {hours_in_day} hours = {hourly_rate:.2f} cc/hr")

calculate_maintenance_fluids()