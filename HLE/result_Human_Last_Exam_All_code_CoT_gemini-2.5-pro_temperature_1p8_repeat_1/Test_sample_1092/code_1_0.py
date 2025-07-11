def calculate_fade_time():
    """
    Calculates the time in years for a Just Noticeable Fade (JNF) to occur
    for a light-sensitive object.
    """
    # Step 1: Define constants based on the problem and conservation standards.
    # Total exposure in lux-hours to cause a JNF for a Bluewool 1 rated material.
    total_lux_hours_for_jnf = 1500000
    
    # Light intensity in lux.
    lux_level = 50
    
    # Assumed daily exposure in hours.
    hours_per_day = 8
    
    # Average number of days in a year to account for leap years.
    days_per_year = 365.25

    # Step 2: Calculate the total daily exposure in lux-hours.
    daily_lux_hours = lux_level * hours_per_day
    
    # Step 3: Calculate the total number of days to reach the JNF threshold.
    total_days_to_fade = total_lux_hours_for_jnf / daily_lux_hours
    
    # Step 4: Convert the total number of days into years.
    years_to_fade = total_days_to_fade / days_per_year

    # Print the explanation and the final equation
    print("To find the time until the next Just Noticeable Fade, we use the following formula:")
    print("Years = (Total Lux-Hours for JNF) / (Lux Level * Hours per Day) / (Days per Year)\n")
    print("Plugging in the values:")
    print(f"Years = {total_lux_hours_for_jnf} / ({lux_level} * {hours_per_day}) / {days_per_year}")
    print(f"\nResult: The next just noticeable fade will occur in approximately {years_to_fade:.2f} years.")

calculate_fade_time()