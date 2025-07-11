import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur
    for a light-sensitive material.
    """
    # 1. Define constants and inputs
    # ISO Bluewool 1 rating corresponds to a Just Noticeable Fade (JNF)
    # at a cumulative exposure of ~20,000 lux-hours.
    lux_hours_for_jnf = 20000

    # Light level the object is exposed to.
    light_level_lux = 50

    # Assuming a standard 8-hour "daily" exposure for a displayed object.
    hours_per_day = 8
    
    # Using 365.25 to account for leap years.
    days_per_year = 365.25

    # 2. Calculate the total annual light exposure
    annual_lux_hours = light_level_lux * hours_per_day * days_per_year
    
    # 3. Calculate the time in years for the fade to occur
    years_to_fade = lux_hours_for_jnf / annual_lux_hours

    # 4. Print the final equation with all the values
    print("The time in years until a Just Noticeable Fade occurs is calculated as:")
    print("Years = (Total Lux-Hours for JNF) / (Lux Level * Hours per Day * Days per Year)")
    print("\nPlugging in the numbers:")
    # The f-string will display each number used in the calculation.
    # We use math.trunc for the annual lux hours to keep the equation clean.
    print(f"Years = {lux_hours_for_jnf} / ({light_level_lux} * {hours_per_day} * {days_per_year})")
    print(f"Years = {lux_hours_for_jnf} / {math.trunc(annual_lux_hours)}")
    
    # Output the final result
    print(f"\nResult: It will take approximately {years_to_fade:.3f} years for the next just noticeable fade to occur.")
    
    # The final answer to be extracted
    global final_answer
    final_answer = f"{years_to_fade:.3f}"

# Run the calculation
calculate_fade_time()
print(f'<<<{final_answer}>>>')