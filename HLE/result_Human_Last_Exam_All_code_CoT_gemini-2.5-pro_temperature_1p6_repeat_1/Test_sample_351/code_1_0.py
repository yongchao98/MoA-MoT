import math

def calculate_maintenance_fluids():
    """
    Calculates maintenance fluid needs for a pediatric burn patient using the BSA method.
    """
    # Patient data
    weight_kg = 25
    tbsa_burn_percent = 45
    bsa_m2 = 0.95

    # Standard BSA method for maintenance fluids: 1500 mL per m^2 per day.
    fluid_per_bsa_per_day = 1500  # mL/m^2/day

    # Step 1: Calculate total daily maintenance fluid.
    daily_fluid_ml = fluid_per_bsa_per_day * bsa_m2

    # Step 2: Convert daily fluid to an hourly rate.
    hourly_fluid_cc_hr = daily_fluid_ml / 24

    # Print the explanation and the calculation steps
    print("This calculation determines the maintenance fluid needs, not the burn resuscitation fluid.")
    print("The Body Surface Area (BSA) method is used (1500 mL/m^2/day).\n")

    print(f"Step 1: Calculate daily maintenance fluid volume based on a BSA of {bsa_m2} m^2.")
    print(f"Equation: {fluid_per_bsa_per_day} mL/m^2 * {bsa_m2} m^2 = {daily_fluid_ml:.1f} mL/day\n")

    print("Step 2: Convert the daily volume to an hourly rate.")
    print(f"Equation: {daily_fluid_ml:.1f} mL / 24 hours = {hourly_fluid_cc_hr:.1f} cc/hr\n")

    print(f"The final calculated maintenance fluid rate is {hourly_fluid_cc_hr:.1f} cc/hr.")

    # Return the final numerical answer
    return round(hourly_fluid_cc_hr, 1)

if __name__ == '__main__':
    final_answer = calculate_maintenance_fluids()
    print(f'<<<{final_answer}>>>')
