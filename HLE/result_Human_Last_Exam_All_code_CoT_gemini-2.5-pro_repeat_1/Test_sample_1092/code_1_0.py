import math

def calculate_fading_time():
    """
    Calculates the time in years for a Just Noticeable Fade (JNF) to occur
    for a highly light-sensitive material under specific light conditions.
    """
    # Standard light dose in lux-hours to cause a JNF for ISO Bluewool 1 rated material.
    lux_hours_for_jnf = 50000

    # Given light level in lux.
    illuminance_lux = 50

    # Damage factor for UV-rich light. Assumed to be twice as damaging as standard UV-filtered light.
    uv_damage_factor = 2

    # Assumed daily exposure in hours (standard for museums/galleries).
    hours_per_day = 8

    # Days in a year.
    days_per_year = 365

    # --- Calculations ---

    # 1. Calculate the effective lux by applying the UV damage factor.
    effective_lux = illuminance_lux * uv_damage_factor

    # 2. Calculate the total hours of exposure needed to reach the JNF threshold.
    total_hours_to_fade = lux_hours_for_jnf / effective_lux

    # 3. Calculate the total number of years to reach the JNF.
    years_to_fade = total_hours_to_fade / (hours_per_day * days_per_year)

    # --- Output ---

    print("This calculation determines the time until a 'Just Noticeable Fade' (JNF) occurs.")
    print(f"A material with ISO Bluewool rating 1 fades at a light dose of {lux_hours_for_jnf} lux-hours.")
    print(f"The effective illuminance is {illuminance_lux} lux * {uv_damage_factor} (UV damage factor) = {effective_lux} lux.")
    print(f"Assuming daily exposure of {hours_per_day} hours for {days_per_year} days a year.")
    print("\nThe final equation is:")
    print(f"Years to Fade = {lux_hours_for_jnf} / ({effective_lux} * {hours_per_day} * {days_per_year})")
    print(f"\nThe next just noticeable fade will occur in approximately {years_to_fade:.2f} years.")

calculate_fading_time()