def calculate_burn_maintenance_fluids():
    """
    Calculates the maintenance fluid needs for a pediatric burn patient based on BSA.
    """
    # Patient data
    bsa_m2 = 0.95  # Body Surface Area in m^2

    # Standard values for calculation
    maintenance_rate_per_m2_per_day = 2000  # in ml/m^2/day
    hours_in_a_day = 24

    # --- Plan explanation ---
    print("To calculate the maintenance fluid needs for a pediatric burn patient, we use a formula based on Body Surface Area (BSA) to account for increased metabolic demands.")
    print(f"The standard rate is {maintenance_rate_per_m2_per_day} ml per square meter of BSA per day.")
    print("\nFirst, we calculate the total 24-hour fluid volume:")
    
    # --- Calculation step 1: Total daily volume ---
    total_fluid_24h = maintenance_rate_per_m2_per_day * bsa_m2
    print(f"Equation: {maintenance_rate_per_m2_per_day} ml/m^2/day * {bsa_m2} m^2 = {total_fluid_24h:.1f} ml/day")

    print("\nNext, we divide the total daily volume by 24 hours to find the hourly rate in cc/hr (1 ml = 1 cc).")

    # --- Calculation step 2: Hourly rate ---
    hourly_rate_cc_hr = total_fluid_24h / hours_in_a_day
    
    print("\nFinal Calculation:")
    print(f"({maintenance_rate_per_m2_per_day} * {bsa_m2}) / {hours_in_a_day} = {hourly_rate_cc_hr:.1f} cc/hr")
    
    # The final numerical answer for the system to read.
    global final_answer
    final_answer = round(hourly_rate_cc_hr, 1)


calculate_burn_maintenance_fluids()
<<<79.2>>>