import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur
    for a highly sensitive material.
    """
    # --- Constants and Assumptions ---

    # ISO Bluewool 1 materials show a "just noticeable fade" (JNF)
    # at a cumulative exposure of ~50,000 lux hours.
    total_dose_for_jnf = 50000  # in lux hours

    # Light intensity the object is exposed to.
    light_intensity = 50  # in lux

    # Assumed daily exposure duration in a gallery or display setting.
    hours_per_day = 8  # in hours

    # Average number of days in a year to account for leap years.
    days_per_year = 365.25

    # --- Calculations ---

    # 1. Calculate the daily light exposure dose.
    daily_dose = light_intensity * hours_per_day

    # 2. Calculate the total number of days to reach the JNF threshold.
    days_to_fade = total_dose_for_jnf / daily_dose

    # 3. Convert the total days to years.
    years_to_fade = days_to_fade / days_per_year

    # --- Output ---
    print("This calculation determines the time to the next 'just noticeable fade' based on standard conservation data.")
    print(f"Assumption: The object is lit for {hours_per_day} hours per day.\n")
    print("The formula is: (Total Dose for Fade / (Light Intensity * Hours per Day)) / Days per Year\n")

    print("Plugging in the numbers for the equation:")
    # Using math.trunc to show the integer part for clarity in the equation printout
    print(f"({total_dose_for_jnf} lux-hours / ({light_intensity} lux * {hours_per_day} hours/day)) / {days_per_year} days/year")
    print(f"= ({total_dose_for_jnf} / {daily_dose}) days / {days_per_year} days/year")
    print(f"= {math.trunc(days_to_fade)} days / {days_per_year} days/year\n")

    print(f"The next just noticeable fade will occur in approximately {years_to_fade:.2f} years.")

    # Returning the final value for the specialized output format
    return years_to_fade

# Execute the function and capture the result for the final output
final_answer = calculate_fade_time()
print(f"<<<{final_answer:.2f}>>>")