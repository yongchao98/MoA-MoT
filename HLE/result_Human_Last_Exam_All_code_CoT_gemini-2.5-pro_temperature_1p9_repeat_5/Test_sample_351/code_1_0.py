def calculate_maintenance_fluid():
    """
    Calculates the maintenance fluid needs for a pediatric patient based on Body Surface Area (BSA).
    """
    # Patient data
    bsa_m2 = 0.95  # Body Surface Area in square meters
    
    # Constants for BSA method
    daily_rate_per_m2 = 1500  # Standard rate in mL/m^2/day
    hours_in_day = 24

    # --- Calculation ---
    # 1. Calculate total daily fluid volume
    total_daily_fluid_ml = daily_rate_per_m2 * bsa_m2
    
    # 2. Calculate the hourly rate
    # 1 mL is equal to 1 cc
    hourly_rate_cc_hr = total_daily_fluid_ml / hours_in_day
    
    # --- Output ---
    # Print the equation and the final result
    print("Maintenance fluid calculation based on the Body Surface Area (BSA) method.")
    print("The formula is: (Daily Rate per m^2 * BSA) / 24 hours")
    print(f"({daily_rate_per_m2} mL/m^2 * {bsa_m2} m^2) / {hours_in_day} hr = {hourly_rate_cc_hr:.1f} cc/hr")

calculate_maintenance_fluid()