import math

def calculate_fade_time():
    """
    Calculates the time in years for a Just Noticeable Fade (JNF) to occur
    for a highly light-sensitive material.
    """
    # 1. Define constants and parameters based on the problem and conservation science standards.
    
    # The light level in lux.
    lux_level = 50
    
    # ISO Bluewool Rating for the material.
    iso_bluewool_rating = 1
    
    # Cumulative exposure in lux-hours to cause one JNF for a Bluewool 1 object.
    # This is a standard value from conservation science.
    total_lux_hours_for_jnf = 15000
    
    # Assume a standard museum display day.
    hours_per_day = 8
    
    # Use 365.25 to account for leap years for better accuracy over time.
    days_per_year = 365.25

    # 2. Perform the calculations step-by-step.
    
    # Calculate the total light exposure (dose) per day.
    daily_lux_hours = lux_level * hours_per_day
    
    # Calculate the total number of days it will take to reach the JNF threshold.
    total_days_to_fade = total_lux_hours_for_jnf / daily_lux_hours
    
    # Convert the total days into years.
    years_to_fade = total_days_to_fade / days_per_year

    # 3. Print the results in a clear, step-by-step explanation.
    
    print(f"Problem: For a material with ISO Bluewool Rating {iso_bluewool_rating}, how long until a 'Just Noticeable Fade' (JNF) occurs at {lux_level} lux?\n")
    print("--- Calculation ---")
    
    print(f"\nStep 1: Calculate the daily light dose.")
    print("This is based on an assumed display time of 8 hours per day.")
    # The final equation requires printing each number
    print(f"Equation: {lux_level} lux * {hours_per_day} hours/day")
    print(f"Result: {daily_lux_hours} lux-hours per day\n")

    print("Step 2: Calculate the number of days to cause one JNF.")
    print(f"A JNF for Bluewool 1 occurs at a total dose of ~{total_lux_hours_for_jnf} lux-hours.")
    # The final equation requires printing each number
    print(f"Equation: {total_lux_hours_for_jnf} lux-hours / {daily_lux_hours} lux-hours/day")
    print(f"Result: {total_days_to_fade:.2f} days\n")

    print("Step 3: Convert the total days to years.")
    # The final equation requires printing each number
    print(f"Equation: {total_days_to_fade:.2f} days / {days_per_year} days/year")
    print(f"Result: {years_to_fade:.2f} years\n")
    
    print("--- Conclusion ---")
    print(f"It will take approximately {years_to_fade:.2f} years for the next just noticeable fade to occur.")
    print("Note: The 'UV-rich' light specified would likely accelerate this process further.")

calculate_fade_time()