def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur
    for a highly sensitive object under specific light conditions.
    """
    # 1. Define constants and variables based on the problem and conservation standards.
    
    # Fading threshold in lux hours for a material with ISO Bluewool Rating 1.
    lux_hours_for_jnf = 20000.0

    # Light exposure conditions from the problem.
    lux_level = 50.0

    # Assumed standard daily duration of light exposure in hours.
    hours_per_day = 8.0

    # A standard damage multiplier for UV-rich light sources.
    uv_damage_factor = 3.0
    
    # Number of days in a year, accounting for leap years.
    days_per_year = 365.25

    # 2. Calculate the effective annual light exposure.
    effective_annual_exposure = lux_level * hours_per_day * uv_damage_factor * days_per_year

    # 3. Calculate the time in years until a Just Noticeable Fade (JNF).
    years_to_fade = lux_hours_for_jnf / effective_annual_exposure

    # 4. Print the explanation, the full equation with values, and the result.
    print("Based on conservation standards, the time to the next just noticeable fade is calculated by dividing the material's fading threshold (in lux hours) by the effective annual exposure.")
    print("\nHere is the equation using the specified values:")
    
    # The final equation string showing all the numbers.
    equation_string = (
        f"Years to Fade = {int(lux_hours_for_jnf)} / "
        f"({int(lux_level)} lux * {int(hours_per_day)} hours/day * "
        f"{int(uv_damage_factor)} UV_factor * {days_per_year} days/year)"
    )
    
    print(equation_string)

    print(f"\nThe object will experience a just noticeable fade in approximately {years_to_fade:.4f} years.")

calculate_fade_time()