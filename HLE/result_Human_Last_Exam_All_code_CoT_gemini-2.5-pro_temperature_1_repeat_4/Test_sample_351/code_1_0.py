def calculate_pediatric_burn_fluid_rate():
    """
    Calculates the total fluid rate for a pediatric burn patient.

    This calculation combines the Parkland formula for burn resuscitation fluid
    with a BSA-based formula for maintenance fluid to determine the total
    fluid requirement over 24 hours, then provides an average hourly rate.
    """
    # Patient data
    weight_kg = 25.0
    tbsa_percent = 45.0
    bsa_m2 = 0.95

    # Constants for the formulas
    parkland_multiplier = 4.0  # mL
    maintenance_multiplier = 1500.0  # mL/m^2
    hours_in_day = 24.0

    # Step 1: Calculate the 24-hour resuscitation fluid volume (Parkland formula)
    resuscitation_fluid_24hr = parkland_multiplier * weight_kg * tbsa_percent
    print("Step 1: Calculate 24-hour resuscitation fluid volume (mL).")
    print(f"Formula: {parkland_multiplier} mL * {weight_kg} kg * {tbsa_percent} %TBSA")
    print(f"Result: {resuscitation_fluid_24hr} mL\n")

    # Step 2: Calculate the 24-hour maintenance fluid volume (BSA method)
    maintenance_fluid_24hr = maintenance_multiplier * bsa_m2
    print("Step 2: Calculate 24-hour maintenance fluid volume (mL).")
    print(f"Formula: {maintenance_multiplier} mL/m^2 * {bsa_m2} m^2")
    print(f"Result: {maintenance_fluid_24hr} mL\n")

    # Step 3: Calculate the total fluid volume needed over 24 hours
    total_fluid_24hr = resuscitation_fluid_24hr + maintenance_fluid_24hr
    print("Step 3: Calculate total 24-hour fluid volume (mL).")
    print(f"Formula: {resuscitation_fluid_24hr} mL + {maintenance_fluid_24hr} mL")
    print(f"Result: {total_fluid_24hr} mL\n")

    # Step 4: Calculate the average hourly rate in cc/hr (1 mL = 1 cc)
    hourly_rate_cc_hr = total_fluid_24hr / hours_in_day
    print("Step 4: Calculate the final fluid rate in cc/hr.")
    print(f"Formula: {total_fluid_24hr} mL / {hours_in_day} hours")
    print(f"Result: {hourly_rate_cc_hr} cc/hr\n")

    print(f"The final calculated fluid rate for the patient is {hourly_rate_cc_hr:.3f} cc/hr.")
    return hourly_rate_cc_hr

# Execute the calculation and store the final answer
final_answer = calculate_pediatric_burn_fluid_rate()
# The final answer is wrapped below as requested
# print(f'<<<{final_answer:.3f}>>>')