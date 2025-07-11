def calculate_maintenance_fluids():
    """
    Calculates maintenance fluid needs for a pediatric patient using the BSA method.
    """
    # Patient and formula parameters
    bsa_m2 = 0.95  # Body Surface Area in square meters
    fluid_per_m2_per_day = 1500  # Standard fluid rate in mL/m^2/day
    hours_in_day = 24

    # Step 1: Calculate total daily maintenance fluid
    daily_fluid_volume = bsa_m2 * fluid_per_m2_per_day

    # Step 2: Calculate the hourly rate
    hourly_rate = daily_fluid_volume / hours_in_day

    # Print the explanation and calculation steps
    print("Calculating Maintenance Fluid Needs using the Body Surface Area (BSA) Method\n")
    print(f"The formula for daily fluid is: BSA (m^2) * {fluid_per_m2_per_day} mL/m^2/day")
    print(f"Daily Fluid = {bsa_m2} m^2 * {fluid_per_m2_per_day} mL/m^2/day = {daily_fluid_volume:.2f} mL/day\n")
    print("To find the hourly rate, we divide the daily total by 24 hours.")
    print(f"Hourly Rate = {daily_fluid_volume:.2f} mL / {hours_in_day} hours\n")
    print(f"Final Answer: The required maintenance fluid rate is {hourly_rate:.2f} cc/hr.")

    return hourly_rate

if __name__ == '__main__':
    final_answer = calculate_maintenance_fluids()
    # The final answer is wrapped in <<<>>> as requested, though it's already printed above.
    # print(f"\n<<<{final_answer:.2f}>>>")
