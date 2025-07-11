import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    for a highly light-sensitive material (ISO Bluewool 1).
    """

    # Constants based on conservation science standards
    # Total light exposure (in lux-hours) to cause a JNF for Bluewool 1.
    # The standard is 0.1 Megalux-hours, which equals 100,000 lux-hours.
    total_exposure_for_jnf = 100000.0

    # Conditions specified in the problem
    lux_level = 50.0

    # A standard assumption for museum and gallery display hours
    hours_per_day = 8.0
    days_per_year = 365.0

    # Step 1: Calculate the annual light exposure in lux-hours per year.
    annual_exposure = lux_level * hours_per_day * days_per_year

    # Step 2: Calculate the time in years to reach the JNF threshold.
    years_to_fade = total_exposure_for_jnf / annual_exposure

    # Step 3: Print the explanation, the equation with all numbers, and the final result.
    print("The time to the next just noticeable fade is calculated as follows:")
    print("We assume a standard museum display time of 8 hours per day.")
    print("\n--- Calculation Breakdown ---")
    print(f"Total Exposure for Fade = {int(total_exposure_for_jnf)} lux-hours (Standard for Bluewool 1)")
    print(f"Annual Exposure = {int(lux_level)} lux * {int(hours_per_day)} hours/day * {int(days_per_year)} days/year = {int(annual_exposure)} lux-hours/year")
    print("\n--- Final Equation ---")
    print(f"Time to Fade (Years) = Total Exposure for Fade / Annual Exposure")
    print(f"Time to Fade (Years) = {int(total_exposure_for_jnf)} / {int(annual_exposure)}")
    print(f"Result: It will take approximately {years_to_fade:.3f} years for the next just noticeable fade to occur.")

# Run the calculation and print the results
calculate_fade_time()

<<<0.685>>>