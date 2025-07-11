def calculate_pediatric_burn_maintenance_fluids():
    """
    Calculates the maintenance fluid needs for a pediatric burn patient.

    The formula used is:
    Total Rate (cc/hr) = Basal Needs (cc/hr) + Evaporative Water Loss (cc/hr)
    - Basal Needs = (1500 mL/m^2/day * BSA) / 24 hr/day
    - Evaporative Water Loss = (25 + %TBSA) * BSA
    """
    # Patient data
    weight_kg = 25
    tbsa_percent = 45
    bsa_m2 = 0.95

    # Formula constants
    basal_rate_per_m2_per_day = 1500  # mL/m^2/day
    ewl_base_factor = 25  # for the EWL formula (25 + %TBSA) * BSA

    # --- Step 1: Calculate Basal Fluid Needs ---
    daily_basal_needs = basal_rate_per_m2_per_day * bsa_m2
    hourly_basal_needs = daily_basal_needs / 24

    # --- Step 2: Calculate Evaporative Water Loss (EWL) ---
    hourly_ewl = (ewl_base_factor + tbsa_percent) * bsa_m2

    # --- Step 3: Calculate Total Hourly Rate ---
    total_hourly_rate = hourly_basal_needs + hourly_ewl

    # --- Output the results ---
    print("Calculating Pediatric Burn Maintenance Fluid Rate:")
    print("-" * 50)
    print(f"Patient Data:")
    print(f"  - Weight: {weight_kg} kg")
    print(f"  - TBSA Burn: {tbsa_percent}%")
    print(f"  - BSA: {bsa_m2} m^2")
    print("\nCalculation Steps:")
    print(f"1. Hourly Basal Needs = ({basal_rate_per_m2_per_day} * {bsa_m2}) / 24 = {hourly_basal_needs:.1f} cc/hr")
    print(f"2. Hourly Evaporative Loss = ({ewl_base_factor} + {tbsa_percent}) * {bsa_m2} = {hourly_ewl:.1f} cc/hr")
    print("\nFinal Calculation:")
    # The request asks to output each number in the final equation.
    print(f"Total Maintenance Rate = {hourly_basal_needs:.1f} cc/hr (Basal) + {hourly_ewl:.1f} cc/hr (Evaporative)")
    print("-" * 50)
    print(f"The final calculated maintenance fluid rate is {total_hourly_rate:.1f} cc/hr.")


if __name__ == "__main__":
    calculate_pediatric_burn_maintenance_fluids()