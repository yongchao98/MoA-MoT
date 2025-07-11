import math

def calculate_fading_time():
    """
    Calculates the time in years for a just noticeable fade to occur
    for a highly light-sensitive object.
    """
    # --- Conservation Science Constants ---
    # Total light exposure in lux-hours to cause a Just Noticeable Fade (JNF)
    # for a material with an ISO Bluewool rating of 1.
    total_lux_hours_for_jnf = 100000.0

    # --- User-Provided Parameters ---
    # Light intensity in lux.
    lux_intensity = 50
    # A standard 8-hour display day is assumed as "daily" exposure.
    hours_per_day = 8
    # Number of days in a year.
    days_per_year = 365

    # --- Calculation ---
    # Step 1: Calculate the total light exposure the object receives per year.
    annual_exposure = lux_intensity * hours_per_day * days_per_year

    # Step 2: Calculate the number of years until a Just Noticeable Fade occurs.
    years_to_fade = total_lux_hours_for_jnf / annual_exposure

    # --- Output ---
    print("This script calculates the time until the next just noticeable fade occurs.")
    print("The calculation is: Years = (Total Lux-Hours for JNF) / (Annual Lux Exposure)")
    print("\nFirst, we calculate the Annual Lux Exposure:")
    print(f"Annual Exposure = {lux_intensity} lux * {hours_per_day} hours/day * {days_per_year} days/year")
    print(f"Annual Exposure = {annual_exposure} lux-hours/year")

    print("\nNext, we use this value to calculate the years to fade for a Bluewool 1 rated material:")
    print("The final equation is:")
    print(f"Years to Fade = {total_lux_hours_for_jnf} / {annual_exposure}")

    print(f"\nResult: The next just noticeable fade will occur in approximately {years_to_fade:.2f} years.")
    
    # Returning the final numeric value for the answer tag.
    return years_to_fade

# Run the calculation
final_answer = calculate_fading_time()