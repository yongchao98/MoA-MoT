def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    on a light-sensitive object based on its Bluewool rating and light exposure.
    """

    # 1. Define constants based on conservation standards and problem statement.
    # jnf_dose_bw1 is the cumulative exposure in lux-hours for a JNF on a
    # material with an ISO Bluewool rating of 1. This is a standard value.
    jnf_dose_bw1 = 50000

    # lux_level is the given light intensity.
    lux_level = 50

    # We assume a standard daily exposure time of 8 hours for 365 days.
    hours_per_day = 8
    days_per_year = 365

    # 2. Calculate the total annual light exposure.
    # Annual Exposure = Lux Level * Hours per Day * Days per Year
    annual_exposure = lux_level * hours_per_day * days_per_year

    # 3. Calculate the number of years until a JNF occurs.
    # Years to Fade = Total JNF Dose / Annual Exposure
    years_to_fade = jnf_dose_bw1 / annual_exposure

    # Print the final equation showing all the numbers involved.
    print(f"The calculation for time to fade (in years) is:")
    print(f"Total JNF Dose / (Lux Level * Daily Hours * Days per Year)")
    print(f"\nEquation with values:")
    print(f"{jnf_dose_bw1} / ({lux_level} * {hours_per_day} * {days_per_year}) = {years_to_fade:.2f} years")

# Execute the function to print the result.
calculate_fade_time()