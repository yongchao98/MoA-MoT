import math

def calculate_fading_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    for a highly light-sensitive object (ISO Bluewool 1).
    """

    # --- Constants based on conservation standards ---
    
    # Total light dose (in lux-hours) to cause one JNF for an
    # ISO Bluewool 1 rated material. This is a standard value.
    total_lux_hours_for_jnf = 200000

    # Light intensity from the user's query.
    lux = 50

    # Assuming a standard exhibition day of 8 hours for "daily" exposure.
    hours_per_day = 8

    # Number of days in a year.
    days_per_year = 365

    # --- Calculation ---
    
    # 1. Calculate the total annual light exposure in lux-hours.
    annual_lux_hours = lux * hours_per_day * days_per_year

    # 2. Calculate the number of years to reach the JNF threshold.
    years_to_fade = total_lux_hours_for_jnf / annual_lux_hours

    # --- Output ---
    
    print("The time to the next just noticeable fade is calculated as:")
    print("Years = (Total Lux-Hours for JNF) / (Annual Lux-Hours)\n")
    
    print("Plugging in the numbers for the final equation:")
    # Print the equation with all the numbers
    print(f"Years = {total_lux_hours_for_jnf} / ({lux} * {hours_per_day} * {days_per_year})")
    
    print(f"\nWhich simplifies to:")
    print(f"Years = {total_lux_hours_for_jnf} / {annual_lux_hours}")

    print(f"\nFinal Result: The next just noticeable fade will occur in approximately {years_to_fade:.2f} years.")

if __name__ == "__main__":
    calculate_fading_time()

<<<1.37>>>