def calculate_maintenance_fluids():
    """
    Calculates maintenance fluid needs for a pediatric patient using the BSA method.
    """
    # Patient information
    bsa_m2 = 0.95  # Body Surface Area in m^2

    # Constants
    # Standard maintenance fluid rate is 1500 ml per square meter of BSA per day.
    maintenance_rate_per_bsa = 1500  # ml/m^2/day
    hours_in_day = 24

    # Step 1: Calculate total daily maintenance fluid
    daily_fluid_ml = maintenance_rate_per_bsa * bsa_m2

    # Step 2: Convert daily fluid needs to an hourly rate
    hourly_fluid_cc = daily_fluid_ml / hours_in_day

    # --- Output of the calculation steps ---
    print("This calculation determines the baseline maintenance fluid needs, not the additional resuscitation fluids required for the burn itself.")
    print("-" * 50)
    print("Step 1: Calculate Daily Maintenance Fluid (ml/day)")
    print(f"Formula: Maintenance Rate (ml/m^2/day) * BSA (m^2)")
    print(f"Calculation: {maintenance_rate_per_bsa} * {bsa_m2} = {daily_fluid_ml:.2f} ml/day")
    print("-" * 50)
    print("Step 2: Convert Daily Fluid to Hourly Rate (cc/hr)")
    print(f"Formula: Daily Fluid (ml) / 24 hours")
    # Note: 1 ml = 1 cc
    print(f"Calculation: {daily_fluid_ml:.2f} / {hours_in_day} = {hourly_fluid_cc:.3f} cc/hr")
    print("-" * 50)
    print(f"The final calculated maintenance fluid rate is {hourly_fluid_cc:.3f} cc/hr.")

calculate_maintenance_fluids()