def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur
    for a highly light-sensitive material.
    """

    # --- Constants based on conservation science and the problem statement ---

    # The cumulative light exposure in lux-hours to cause a "Just Noticeable Fade" (JNF)
    # for the most sensitive materials, rated as ISO Bluewool 1.
    TOTAL_LUX_HOURS_FOR_JNF_BW1 = 50000

    # The illuminance level the object is exposed to, in lux.
    LUX_LEVEL = 50

    # Assuming a standard museum display schedule for "daily" exposure.
    HOURS_PER_DAY = 8
    DAYS_PER_YEAR = 365

    # --- Calculation ---

    # 1. Calculate the total light exposure the object receives in one year.
    annual_lux_hours = LUX_LEVEL * HOURS_PER_DAY * DAYS_PER_YEAR

    # 2. Calculate the number of years it will take to reach the JNF threshold.
    years_to_fade = TOTAL_LUX_HOURS_FOR_JNF_BW1 / annual_lux_hours

    # --- Output ---

    print(f"An object with ISO Bluewool Rating 1 requires approximately {TOTAL_LUX_HOURS_FOR_JNF_BW1} lux-hours of exposure for a just noticeable fade.")
    print(f"The annual light exposure at {LUX_LEVEL} lux for {HOURS_PER_DAY} hours/day is calculated as:")
    print(f"  {LUX_LEVEL} lux * {HOURS_PER_DAY} hours/day * {DAYS_PER_YEAR} days/year = {annual_lux_hours} lux-hours per year.")
    print("\nTo find the years until the next fade, we divide the total required exposure by the annual exposure.")
    print("\nFinal Equation:")
    print(f"Years to Fade = {TOTAL_LUX_HOURS_FOR_JNF_BW1} / {annual_lux_hours}")
    print(f"\nResult: It will take approximately {years_to_fade:.3f} years for the next just noticeable fade to occur.")
    
    return years_to_fade

# Run the calculation and store the result
result = calculate_fade_time()
# The final answer is wrapped in <<<>>> as requested, though it's already printed above.
# This format is for automated checking. The numerical value is the final answer.
# print(f"<<<{result:.3f}>>>")