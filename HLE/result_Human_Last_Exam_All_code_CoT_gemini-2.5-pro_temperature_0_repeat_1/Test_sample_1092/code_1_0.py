import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur on a
    highly light-sensitive material under specific light conditions.
    """
    # Step 1: Define constants and inputs based on the problem description.
    # The total cumulative exposure that causes a Just Noticeable Fade (JNF)
    # for a material with ISO Bluewool Rating 1.
    total_lux_hours_for_jnf = 500000

    # The intensity of the light source in lux.
    lux_intensity = 50

    # A damage multiplier for UV-rich light. Standard museum light is filtered,
    # so UV-rich light is more damaging. A factor of 2 is a common estimate.
    uv_damage_factor = 2

    # Assumed hours of light exposure per day (e.g., museum opening hours).
    hours_per_day = 8

    # Average number of days in a year, accounting for leap years.
    days_per_year = 365.25

    # Step 2: Calculate the total effective light exposure per year.
    # This is the effective lux (intensity * UV factor) multiplied by the
    # total hours of exposure in a year.
    effective_annual_exposure = lux_intensity * uv_damage_factor * hours_per_day * days_per_year

    # Step 3: Calculate the number of years until a JNF occurs.
    # This is the total exposure threshold divided by the annual exposure rate.
    years_to_jnf = total_lux_hours_for_jnf / effective_annual_exposure

    # Step 4: Print the explanation and the final equation with all numbers.
    print("The time to the next just noticeable fade is calculated by dividing the total required exposure (in lux-hours) by the annual rate of exposure.")
    print("\nFinal Equation:")
    print(f"Years = (Total Lux-Hours for JNF) / (Lux * UV Factor * Hours per Day * Days per Year)")
    print("\nCalculation:")
    # The f-string below substitutes the variable names with their numerical values to show the exact calculation.
    print(f"Years = {total_lux_hours_for_jnf} / ({lux_intensity} * {uv_damage_factor} * {hours_per_day} * {days_per_year})")
    print(f"Years = {total_lux_hours_for_jnf} / {effective_annual_exposure}")
    print(f"Years = {years_to_jnf:.2f}")


calculate_fade_time()
<<<1.71>>>