import math

def calculate_maintenance_fluids():
    """
    Calculates pediatric maintenance fluid needs using the Body Surface Area (BSA) method.
    """
    # Patient and formula data
    bsa_m2 = 0.95
    daily_rate_per_m2 = 1500  # cc/m^2/day
    hours_per_day = 24

    # --- Calculation ---
    
    # Step 1: Calculate the total daily maintenance fluid requirement.
    daily_fluid_cc = daily_rate_per_m2 * bsa_m2

    # Step 2: Convert the daily requirement to an hourly rate.
    hourly_fluid_rate = daily_fluid_cc / hours_per_day

    # --- Output the results step-by-step ---
    
    print("Calculating maintenance fluid needs using the Body Surface Area (BSA) method.")
    print("Formula: Rate (cc/hr) = (1500 cc/m^2/day * BSA) / 24 hr/day\n")

    print(f"Patient's Body Surface Area (BSA) = {bsa_m2} m^2")
    
    print("\n--- Step 1: Calculate Total Daily Fluid Volume ---")
    # Using math.trunc() for the print output to keep it clean like the problem asks for integers where appropriate.
    # The actual calculation will use the float.
    print(f"Equation: Daily Volume = {daily_rate_per_m2} cc/m^2/day * {bsa_m2} m^2")
    print(f"Result: {daily_fluid_cc} cc/day")
    
    print("\n--- Step 2: Convert Daily Volume to Hourly Rate ---")
    print(f"Equation: Hourly Rate = {daily_fluid_cc} cc/day / {hours_per_day} hours")
    # round the final answer for readability.
    print(f"Result: {hourly_fluid_rate:.2f} cc/hr")

    print(f"\nThe final calculated maintenance fluid rate is {hourly_fluid_rate:.2f} cc/hr.")
    # This last part is just for the final answer extraction format.
    # print(f"<<<{hourly_fluid_rate:.2f}>>>")

calculate_maintenance_fluids()