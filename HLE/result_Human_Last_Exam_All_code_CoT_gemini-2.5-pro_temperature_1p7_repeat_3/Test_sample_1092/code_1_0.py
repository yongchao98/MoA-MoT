def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    on an object with an ISO Bluewool Rating of 1 under specific lighting conditions.
    """

    # Step 1: Define the known constants based on conservation standards.
    # Total light exposure (in lux-hours) to cause one JNF for a Bluewool 1 material.
    total_exposure_for_jnf_bw1 = 7500
    # Light intensity the object is exposed to.
    light_intensity_lux = 50
    # Assuming a standard 8-hour exposure for a 'daily' display.
    hours_per_day = 8
    # Days in a year.
    days_per_year = 365

    # Step 2: Calculate the total annual light exposure.
    annual_exposure_lux_hours = light_intensity_lux * hours_per_day * days_per_year

    # Step 3: Calculate the time in years to reach the JNF threshold.
    years_to_jnf = total_exposure_for_jnf_bw1 / annual_exposure_lux_hours

    # Print the explanation and the steps of the calculation.
    print("To find the time until a 'just noticeable fade' occurs, we calculate the cumulative light exposure.")
    print(f"\n1. Standard for Just Noticeable Fade (JNF) for ISO Bluewool 1: {total_exposure_for_jnf_bw1} lux-hours.")
    
    print("\n2. Calculate the total annual light exposure:")
    print(f"   Equation: Light Intensity * Hours per Day * Days per Year")
    print(f"   Calculation: {light_intensity_lux} lux * {hours_per_day} hours/day * {days_per_year} days/year = {annual_exposure_lux_hours} lux-hours/year")

    print("\n3. Calculate the years until a JNF occurs:")
    print(f"   Equation: Total Exposure for JNF / Annual Exposure")
    print(f"   Final Calculation: {total_exposure_for_jnf_bw1} / {annual_exposure_lux_hours} = {years_to_jnf:.4f} years")

    # The result shows that a just noticeable fade will occur in approximately 0.0514 years,
    # which is equivalent to about 19 days.

calculate_fade_time()