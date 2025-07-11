def calculate_pediatric_burn_fluids():
    """
    Calculates maintenance fluid needs for a pediatric burn patient using the Galveston Formula.
    """
    # Patient data
    weight_kg = 25
    tbsa_percent = 45
    bsa_m2 = 0.95

    # --- Calculation Plan ---

    # 1. Convert TBSA percent to a decimal for calculation
    tbsa_decimal = tbsa_percent / 100

    # 2. Calculate the surface area of the burn in m^2
    burned_area_m2 = bsa_m2 * tbsa_decimal

    # 3. Calculate the resuscitation fluid component (5000 mL per m^2 of burned area)
    fluid_for_burn = 5000 * burned_area_m2

    # 4. Calculate the maintenance fluid component (2000 mL per m^2 of total BSA)
    maintenance_fluid = 2000 * bsa_m2

    # 5. Calculate the total fluid needed over 24 hours
    total_fluid_24hr = fluid_for_burn + maintenance_fluid

    # 6. Calculate the average hourly rate in cc/hr (1 mL = 1 cc)
    hourly_rate = total_fluid_24hr / 24

    # --- Output Results ---
    print("Calculating Pediatric Burn Fluid Needs using the Galveston Formula")
    print("Formula: ((5000 mL * Burned Area in m^2) + (2000 mL * Total BSA in m^2)) / 24 hrs")
    print("-" * 30)
    print(f"Patient Weight: {weight_kg} kg")
    print(f"Total Body Surface Area (BSA): {bsa_m2} m^2")
    print(f"Burn Percentage (TBSA): {tbsa_percent}%")
    print("-" * 30)
    print("Step 1: Calculate Total 24-Hour Fluid Volume")
    print(f"  Fluid for Burn: 5000 mL * ({bsa_m2} m^2 * {tbsa_decimal}) = {fluid_for_burn:.1f} mL")
    print(f"  Maintenance Fluid: 2000 mL * {bsa_m2} m^2 = {maintenance_fluid:.1f} mL")
    print(f"  Total Fluid in 24 hours = {fluid_for_burn:.1f} mL + {maintenance_fluid:.1f} mL = {total_fluid_24hr:.1f} mL")
    print("-" * 30)
    print("Step 2: Calculate Hourly Infusion Rate")
    # This line fulfills the requirement to output each number in the final equation
    print(f"  Final Equation: ({fluid_for_burn:.1f} + {maintenance_fluid:.1f}) / 24 = {hourly_rate:.2f} cc/hr")
    print("-" * 30)
    print(f"The required maintenance fluid rate is: {hourly_rate:.2f} cc/hr")

    # The final numerical answer for the platform
    return hourly_rate

# Run the calculation and print the final result
final_answer = calculate_pediatric_burn_fluids()
# <<<168.23>>>