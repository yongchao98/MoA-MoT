import math

def calculate_maintenance_fluids():
    """
    Calculates pediatric maintenance fluid needs based on Body Surface Area (BSA).
    
    This function uses the standard formula of 1500 mL/m^2/day.
    """
    
    # Patient data
    bsa_m2 = 0.95  # Body Surface Area in square meters
    
    # Constants
    daily_rate_per_m2 = 1500  # Standard maintenance fluid rate in mL/m^2/day
    hours_in_day = 24
    
    print("Calculating Pediatric Maintenance Fluid Needs using the BSA Method")
    print("--------------------------------------------------------------")
    print(f"Patient's Body Surface Area (BSA): {bsa_m2} m^2")
    
    # Step 1: Calculate total daily fluid requirement
    total_daily_fluid_ml = daily_rate_per_m2 * bsa_m2
    
    print("\nStep 1: Calculate total daily fluid volume in mL.")
    print(f"Equation: Daily Volume = Daily Rate per m^2 * BSA")
    print(f"Calculation: {daily_rate_per_m2} mL/m^2 * {bsa_m2} m^2 = {total_daily_fluid_ml:.2f} mL/day")

    # Step 2: Convert daily fluid requirement to hourly rate
    hourly_rate_cc_hr = total_daily_fluid_ml / hours_in_day
    
    print("\nStep 2: Convert daily volume to an hourly rate in cc/hr (1 mL = 1 cc).")
    print(f"Equation: Hourly Rate = Daily Volume / Hours in a Day")
    print(f"Calculation: {total_daily_fluid_ml:.2f} mL / {hours_in_day} hours = {hourly_rate_cc_hr:.2f} cc/hr")

    print("\nNote: This calculation is for maintenance fluids only and does not include the additional resuscitation fluids required for the burn itself.")
    
    # The final answer to be extracted
    return hourly_rate_cc_hr

if __name__ == '__main__':
    final_answer = calculate_maintenance_fluids()
    # The prompt asks for the answer to be returned at the end in a specific format.
    # To follow that, we print the final numeric answer separately.
    # The instructions also said not to ask the user to copy/paste the result.
    # So the value will just be part of the final print out, and also available for the <<<>>> tag.
    print(f"\nThe final calculated maintenance fluid rate is {final_answer:.2f} cc/hr.")
    # The problem description requests the final answer in a specific format.
    # For example, <<<9.8>>>
    # Let's format the answer as requested, rounded to two decimal places.
    print(f"\n<<<>>>")
    print(f"<<<{final_answer:.2f}>>>")
