import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    on a light-sensitive object based on its ISO Bluewool rating and light exposure.
    """
    # Step 1: Define constants and inputs based on the problem and conservation standards.
    # JNF threshold for ISO Bluewool 1 is ~50,000 lux-hours.
    total_lux_hours_for_jnf = 50000
    
    # Light exposure level in lux.
    lux_level = 50
    
    # Assumed daily exposure time for a display object.
    hours_per_day = 8
    
    # Days in a year, accounting for leap years.
    days_per_year = 365.25

    # Step 2: Calculate the total annual light exposure in lux-hours.
    annual_exposure = lux_level * hours_per_day * days_per_year

    # Step 3: Calculate the time in years for the JNF to occur.
    years_to_fade = total_lux_hours_for_jnf / annual_exposure

    # Print the explanation and the final equation with all numbers.
    print("The time in years for a just noticeable fade is calculated by dividing the material's total JNF light dose limit by its annual light exposure.\n")
    print("Formula: Years = Total_JNF_Lux_Hours / (Lux_Level * Hours_per_Day * Days_per_Year)\n")
    print("Calculation with the given values:")
    print(f"Years = {total_lux_hours_for_jnf} / ({lux_level} * {hours_per_day} * {days_per_year})")
    print(f"Years = {total_lux_hours_for_jnf} / {annual_exposure}")
    print(f"Years ≈ {years_to_fade:.3f}\n")
    print(f"Therefore, the next just noticeable fade will occur in approximately {years_to_fade:.3f} years.")

calculate_fade_time()
<<<0.342>>>