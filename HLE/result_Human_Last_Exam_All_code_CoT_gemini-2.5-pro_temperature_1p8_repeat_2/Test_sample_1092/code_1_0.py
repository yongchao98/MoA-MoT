import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    for a highly light-sensitive object.
    """
    # --- Constants and Assumptions ---

    # The cumulative light exposure in lux-hours that causes a JNF for an
    # object with an ISO Bluewool rating of 1.
    jnf_threshold_lux_hours = 50000

    # The light level the object is exposed to, in lux.
    light_level_lux = 50

    # The assumed number of hours the object is exposed to light per day.
    # A standard 8-hour museum day is a common assumption.
    hours_per_day = 8

    # The average number of days in a year to account for leap years.
    days_per_year = 365.25

    # --- Calculation ---

    # Step 1: Calculate the total number of hours of exposure needed to reach the JNF threshold.
    total_hours_needed = jnf_threshold_lux_hours / light_level_lux

    # Step 2: Calculate how many days of exposure this corresponds to.
    total_days_needed = total_hours_needed / hours_per_day

    # Step 3: Convert the total days into years.
    total_years = total_days_needed / days_per_year

    # --- Output ---
    
    print("This script calculates the time until a Just Noticeable Fade (JNF) occurs based on standard conservation values.")
    print("\nKey values used in the calculation:")
    print(f"- JNF Threshold for ISO Bluewool 1: {jnf_threshold_lux_hours} lux-hours")
    print(f"- Light Intensity: {light_level_lux} lux")
    print(f"- Assumed Daily Exposure: {hours_per_day} hours/day")
    print(f"- Days per Year: {days_per_year}\n")
    
    print("The final equation is structured as follows:")
    print("Years = (JNF Threshold / Light Intensity) / Hours per Day / Days per Year\n")

    print("Plugging in the numbers:")
    # This line prints the final equation with each number, as requested.
    print(f"Years = ({jnf_threshold_lux_hours} / {light_level_lux}) / {hours_per_day} / {days_per_year}")

    print(f"\nResult: {total_years:.2f} years")

# Execute the function to print the result.
calculate_fade_time()