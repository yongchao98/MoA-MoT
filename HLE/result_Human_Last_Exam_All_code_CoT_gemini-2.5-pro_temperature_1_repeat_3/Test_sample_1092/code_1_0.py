import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    on a light-sensitive object based on museum conservation standards.
    """
    # Standard cumulative exposure in lux-hours to cause one JNF for ISO Bluewool 1 material.
    total_lux_hours_for_jnf = 50000

    # Given exposure conditions.
    lux_level = 50

    # Assuming a standard museum display day of 8 hours.
    hours_per_day = 8
    
    # Using 365.25 to account for leap years for better accuracy.
    days_per_year = 365.25

    # Calculate the annual light exposure.
    annual_exposure = lux_level * hours_per_day * days_per_year

    # Calculate the number of years until one JNF occurs.
    years_to_jnf = total_lux_hours_for_jnf / annual_exposure

    print("The time for a 'just noticeable fade' is calculated by dividing the total required lux-hours for fading by the annual light exposure.")
    print(f"Standard for ISO Bluewool 1: {total_lux_hours_for_jnf} lux-hours for one fade.")
    print(f"Exposure conditions: {lux_level} lux for {hours_per_day} hours a day.")
    print("\nFinal Equation:")
    print(f"Years to Fade = {total_lux_hours_for_jnf} / ({lux_level} * {hours_per_day} * {days_per_year})")
    print(f"\nResult: The object will show a just noticeable fade in approximately {years_to_jnf:.3f} years.")

    # Return the raw number for the final answer block.
    return years_to_jnf

# Run the calculation and store the result.
result = calculate_fade_time()
# The final answer will be formatted and appended at the end by the system.
# print(f"\n<<<{result:.3f}>>>") # for final answer formatting.