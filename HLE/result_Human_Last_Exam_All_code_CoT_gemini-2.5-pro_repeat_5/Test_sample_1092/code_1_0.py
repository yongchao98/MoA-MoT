import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    for a highly sensitive material under specific lighting conditions.
    """
    # Step 1: Define the constants and given values.

    # Total light exposure (in lux-hours) to cause a JNF for ISO Bluewool 1 materials.
    # This is a standard value in conservation science.
    total_exposure_for_jnf = 50000  # lux-hours

    # Given light level from the problem.
    lux_level = 50  # lux

    # Assuming a standard 8-hour day for "daily" exposure.
    hours_per_day = 8  # hours

    # Number of days in a year.
    days_per_year = 365

    # A conservative damage factor for "UV-rich" light.
    # This means the light is considered twice as damaging as UV-filtered light.
    uv_damage_factor = 2

    # Step 2: Calculate the effective annual exposure.
    annual_effective_exposure = lux_level * hours_per_day * days_per_year * uv_damage_factor

    # Step 3: Calculate the time in years to reach the JNF threshold.
    years_to_jnf = total_exposure_for_jnf / annual_effective_exposure

    # Step 4: Print the explanation and the final result.
    print("The time to the next just noticeable fade is calculated by dividing the total allowed exposure by the effective annual exposure.")
    print("\n--- Calculation Breakdown ---")
    print(f"Total Allowed Exposure for Bluewool 1: {total_exposure_for_jnf} lux-hours")
    print(f"Light Level: {lux_level} lux")
    print(f"Exposure Duration: {hours_per_day} hours/day for {days_per_year} days/year")
    print(f"UV Damage Factor: {uv_damage_factor}")

    print("\n--- Final Equation ---")
    final_equation = f"Years to Fade = {total_exposure_for_jnf} / ({lux_level} * {hours_per_day} * {days_per_year} * {uv_damage_factor})"
    print(final_equation)

    final_calculation = f"Years to Fade = {total_exposure_for_jnf} / {annual_effective_exposure}"
    print(final_calculation)

    print(f"\nResult: The next just noticeable fade will occur in approximately {years_to_jnf:.2f} years.")

calculate_fade_time()
<<<0.17>>>