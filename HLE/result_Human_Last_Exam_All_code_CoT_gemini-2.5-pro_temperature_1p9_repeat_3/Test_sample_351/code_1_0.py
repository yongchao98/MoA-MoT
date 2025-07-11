import math

def calculate_pediatric_burn_maintenance_fluids():
    """
    Calculates the maintenance fluid needs for a pediatric burn patient using the BSA method.
    """
    # Patient Data
    bsa = 0.95  # Body Surface Area in m^2
    maintenance_factor = 1500  # Standard rate in mL/m^2 per day

    print("This program calculates maintenance fluid needs for a pediatric burn patient using the Body Surface Area (BSA) method.")
    print("The formula is: Daily Fluids (mL) = 1500 mL/m^2 * Patient's BSA (m^2)\n")

    # Step 1: Calculate the total daily maintenance fluid requirement.
    daily_fluid_volume = maintenance_factor * bsa
    print("Step 1: Calculate the total daily fluid volume.")
    print(f"The equation is: {maintenance_factor} mL/m^2 * {bsa} m^2")
    print(f"Result: {daily_fluid_volume:.2f} mL/day\n")

    # Step 2: Convert the daily volume to an hourly rate.
    hours_in_day = 24
    hourly_fluid_rate = daily_fluid_volume / hours_in_day
    print("Step 2: Convert the daily volume to an hourly rate in cc/hr (1 mL = 1 cc).")
    print(f"The equation is: {daily_fluid_volume:.2f} mL / {hours_in_day} hours")
    # Using math.ceil to round up to the nearest whole number as is common practice, but for precision we will show decimals.
    print(f"Result: {hourly_fluid_rate:.2f} cc/hr\n")

    print(f"The final calculated maintenance fluid rate is {hourly_fluid_rate:.2f} cc/hr.")

    return hourly_fluid_rate

# Run the calculation and store the final answer
final_answer = calculate_pediatric_burn_maintenance_fluids()
# The 'final answer' is programmatically generated but we will format it manually as requested.
# For the purpose of the final return value requested by the user prompt
final_answer_for_submission = 59.38
# The following line is for the final answer extraction and is not printed to the user console.
# print(f'<<<{final_answer_for_submission}>>>')