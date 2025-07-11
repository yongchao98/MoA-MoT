def calculate_pediatric_burn_fluids():
    """
    Calculates the initial fluid rate for a pediatric burn patient using the Galveston formula.
    """
    # Patient data
    weight_kg = 25
    tbsa_percent = 45
    bsa_m2 = 0.95

    # Galveston formula constants
    fluid_per_m2_burn = 5000  # mL/m^2 of burn
    fluid_per_m2_maintenance = 2000  # mL/m^2 of total BSA

    print("Calculating the total fluid needs for a pediatric burn patient using the Galveston formula.")
    print("This formula combines resuscitation and maintenance fluids.")
    print("-" * 30)
    print(f"Patient's Total Body Surface Area (BSA): {bsa_m2} m^2")
    print(f"Total Body Surface Area (TBSA) Burn: {tbsa_percent}%")
    print("-" * 30)

    # 1. Calculate the burned surface area in m^2
    burned_bsa_m2 = bsa_m2 * (tbsa_percent / 100.0)
    print(f"Step 1: Calculate burned surface area")
    print(f"Equation: {bsa_m2} m^2 * ({tbsa_percent} / 100) = {burned_bsa_m2:.3f} m^2")
    print()

    # 2. Calculate the fluid for the burn area over 24 hours
    resuscitation_fluid_ml = fluid_per_m2_burn * burned_bsa_m2
    print(f"Step 2: Calculate resuscitation fluid volume for 24 hours")
    print(f"Equation: {fluid_per_m2_burn} mL/m^2 * {burned_bsa_m2:.3f} m^2 = {resuscitation_fluid_ml:.1f} mL")
    print()

    # 3. Calculate the maintenance fluid over 24 hours
    maintenance_fluid_ml = fluid_per_m2_maintenance * bsa_m2
    print(f"Step 3: Calculate maintenance fluid volume for 24 hours")
    print(f"Equation: {fluid_per_m2_maintenance} mL/m^2 * {bsa_m2} m^2 = {maintenance_fluid_ml:.1f} mL")
    print()

    # 4. Calculate total fluid for 24 hours
    total_fluid_24hr_ml = resuscitation_fluid_ml + maintenance_fluid_ml
    print(f"Step 4: Calculate total fluid volume for 24 hours")
    print(f"Equation: {resuscitation_fluid_ml:.1f} mL + {maintenance_fluid_ml:.1f} mL = {total_fluid_24hr_ml:.1f} mL")
    print()

    # 5. Calculate the volume to be infused in the first 8 hours (half of the 24hr total)
    fluid_first_8hr_ml = total_fluid_24hr_ml / 2.0
    print(f"Step 5: Calculate fluid volume for the first 8 hours (half of total)")
    print(f"Equation: {total_fluid_24hr_ml:.1f} mL / 2 = {fluid_first_8hr_ml:.1f} mL")
    print()

    # 6. Calculate the hourly rate for the first 8 hours
    rate_cc_per_hr = fluid_first_8hr_ml / 8.0
    print(f"Step 6: Calculate the hourly infusion rate (cc/hr) for the first 8 hours")
    print(f"Final Equation: {fluid_first_8hr_ml:.1f} mL / 8 hours = {rate_cc_per_hr:.2f} cc/hr")
    print("-" * 30)

    # Store the final numerical answer for the platform
    global final_answer
    final_answer = f"<<<{rate_cc_per_hr:.2f}>>>"

# Run the calculation
calculate_pediatric_burn_fluids()