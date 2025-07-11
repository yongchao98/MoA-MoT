def calculate_fade_time():
    """
    Calculates the time in years for a Just Noticeable Fade (JNF) to occur
    for a highly light-sensitive object.
    """
    # Standard total exposure in lux-hours to cause one JNF for ISO Bluewool 1
    jnf_threshold = 50000

    # Given exposure conditions
    lux = 50
    hours_per_day = 8  # Assuming a standard 8-hour day for display
    days_per_year = 365

    # --- Calculations ---

    # 1. Calculate the total light exposure per year
    annual_exposure = lux * hours_per_day * days_per_year

    # 2. Calculate the number of years to reach the JNF threshold
    years_to_fade = jnf_threshold / annual_exposure

    # --- Output the results ---
    print("This calculation determines the time until one 'Just Noticeable Fade' (JNF) occurs.")
    print(f"JNF Threshold for ISO Bluewool 1: {jnf_threshold} lux-hours")
    print(f"Light Intensity: {lux} lux")
    print(f"Assumed Daily Exposure: {hours_per_day} hours")
    print("\n--- Equation ---")
    print(f"Years to Fade = (JNF Threshold) / (Lux * Hours per Day * Days per Year)")
    print(f"Years to Fade = {jnf_threshold} / ({lux} * {hours_per_day} * {days_per_year})")
    print(f"\nThe next just noticeable fade will occur in approximately {years_to_fade:.2f} years.")

calculate_fade_time()