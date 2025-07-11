def calculate_fade_time():
    """
    Calculates the time in years for one Just Noticeable Fade (JNF) to occur
    for a highly light-sensitive object.
    """
    # According to conservation standards, for a material with ISO Bluewool Rating 1,
    # one Just Noticeable Fade occurs after approximately 200 hours of exposure at 50 lux.
    hours_for_one_fade = 200

    # The object is exposed daily, which we interpret as continuous exposure.
    hours_per_day = 24
    days_per_year = 365

    # Calculate the total number of exposure hours in one year.
    total_annual_hours = hours_per_day * days_per_year

    # Calculate the time in years for one fade to occur.
    # This is the hours needed for one fade divided by the annual exposure hours.
    years_per_fade = hours_for_one_fade / total_annual_hours

    # Print the full equation and the result.
    print("Equation: Time in Years = (Hours for One Fade) / (Hours per Day * Days per Year)")
    print(f"Calculation: {hours_for_one_fade} / ({hours_per_day} * {days_per_year}) = {years_per_fade:.4f} years")

calculate_fade_time()