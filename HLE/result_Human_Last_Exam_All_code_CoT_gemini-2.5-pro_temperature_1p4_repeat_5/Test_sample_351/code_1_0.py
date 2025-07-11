def calculate_maintenance_fluids():
    """
    Calculates pediatric maintenance fluid needs based on Body Surface Area (BSA).
    """
    # Patient data
    bsa = 0.95  # Body Surface Area in m^2
    weight_kg = 25 # not used in BSA method
    tbsa_burn_percent = 45 # not used for maintenance fluid calculation

    # Standard formula values
    daily_fluid_rate_per_m2 = 1500  # Standard rate in mL/m^2/day
    hours_in_day = 24

    # --- Calculation ---
    # 1. Calculate total daily maintenance fluid
    total_daily_fluid = daily_fluid_rate_per_m2 * bsa

    # 2. Calculate hourly maintenance fluid rate
    hourly_fluid_rate = total_daily_fluid / hours_in_day

    # --- Output ---
    print("Calculating pediatric maintenance fluid needs using the Body Surface Area (BSA) method.")
    print("This method calculates the baseline daily fluid requirement.")
    print(f"\nPatient BSA: {bsa} m^2")

    print("\nThe calculation for the hourly rate is (Daily Rate per m^2 * BSA) / 24 hours.")
    print("\nFinal Equation:")
    print(f"({daily_fluid_rate_per_m2} * {bsa}) / {hours_in_day}")

    print(f"\nThe required maintenance fluid rate is: {hourly_fluid_rate:.2f} cc/hr")

if __name__ == "__main__":
    calculate_maintenance_fluids()