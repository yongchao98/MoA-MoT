import math

def calculate_maintenance_fluids():
    """
    Calculates maintenance fluid needs for a pediatric patient using the BSA method.
    """
    # Patient data
    bsa_m2 = 0.95  # Body Surface Area in square meters

    # Constants
    fluid_per_m2_per_day = 1500  # Standard maintenance fluid rate in mL/m^2/day
    hours_per_day = 24

    # Calculate total daily maintenance fluid
    daily_fluid_volume = fluid_per_m2_per_day * bsa_m2

    # Calculate hourly maintenance fluid rate
    hourly_rate = daily_fluid_volume / hours_per_day

    # Print the step-by-step calculation
    print("Calculating maintenance fluid rate using the Body Surface Area (BSA) method.")
    print("\nStep 1: Calculate total daily fluid volume.")
    print(f"Equation: {fluid_per_m2_per_day} mL/m^2/day * {bsa_m2} m^2 = {daily_fluid_volume:.2f} mL/day")

    print("\nStep 2: Calculate the hourly rate.")
    print(f"Equation: {daily_fluid_volume:.2f} mL/day / {hours_per_day} hours/day = {hourly_rate:.2f} mL/hr")

    print(f"\nThe required maintenance fluid rate is {hourly_rate:.2f} cc/hr.")

if __name__ == "__main__":
    calculate_maintenance_fluids()
    # The final numerical answer for the submission system.
    final_answer = (1500 * 0.95) / 24
    print(f'<<<{final_answer:.2f}>>>')
