def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur
    for a light-sensitive object under specific lighting conditions.
    """

    # --- Constants and Given Conditions ---

    # According to conservation science, the total cumulative light exposure
    # to cause a "Just Noticeable Fade" (JNF) for a material with an
    # ISO Bluewool Rating of 1 is approximately 500,000 lux-hours.
    total_exposure_for_jnf = 500000

    # The light level the object is exposed to.
    light_level_lux = 50

    # We assume a standard museum exhibition day of 8 hours.
    hours_per_day = 8

    # Days in a year.
    days_per_year = 365

    # --- Calculation ---

    # 1. Calculate the total exposure per year.
    annual_exposure = light_level_lux * hours_per_day * days_per_year

    # 2. Calculate the number of years to reach the JNF threshold.
    years_to_fade = total_exposure_for_jnf / annual_exposure

    # --- Output ---

    print("To find the time in years until the next just noticeable fade, we divide the total required exposure for a fade by the annual exposure.")
    print("\nCalculation breakdown:")
    print(f"Total exposure for a fade (Bluewool 1): {total_exposure_for_jnf} lux-hours")
    print(f"Annual exposure = {light_level_lux} lux * {hours_per_day} hours/day * {days_per_year} days/year = {annual_exposure} lux-hours/year")
    
    print("\nFinal Equation:")
    print(f"{total_exposure_for_jnf} / ({light_level_lux} * {hours_per_day} * {days_per_year}) = {years_to_fade}")
    print(f"\nThis means the next just noticeable fade will occur in approximately {years_to_fade:.1f} years.")

calculate_fade_time()
<<<3.4>>>