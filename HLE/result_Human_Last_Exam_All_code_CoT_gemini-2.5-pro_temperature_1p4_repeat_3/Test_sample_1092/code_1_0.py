import math

def calculate_fading_time():
    """
    Calculates the time in years for a Just Noticeable Fade (JNF) to occur
    for a light-sensitive object under specific lighting conditions.
    """

    # --- Define Parameters based on the problem and conservation science ---

    # Total light exposure (dose) in lux-hours to cause one JNF for a
    # material with ISO Bluewool rating 1.
    jnf_dose_bw1 = 50000

    # Daily light exposure conditions from the user's query.
    lux_level = 50

    # Assuming a 10-hour display day for "daily" exposure.
    hours_per_day = 10
    
    # Days in a year.
    days_per_year = 365

    # UV damage factor. "UV-rich" light (e.g., unfiltered daylight/fluorescent)
    # is more damaging. A factor of 3 is a standard conservative estimate.
    uv_damage_factor = 3

    # --- Perform Calculation ---

    # 1. Calculate the effective annual exposure in lux-hours.
    # This is the lux level multiplied by the hours of exposure and the UV factor.
    effective_annual_exposure = lux_level * hours_per_day * days_per_year * uv_damage_factor

    # 2. Calculate the time in years to reach the JNF threshold.
    years_to_fade = jnf_dose_bw1 / effective_annual_exposure

    # --- Output the Results ---
    
    print("To find the time to the next Just Noticeable Fade, we use the following equation:")
    print("Years to Fade = Total Fading Dose / (Lux * Hours/Day * Days/Year * UV Factor)\n")
    print("The final calculation with the given values is:")
    
    # Print the equation with all the numbers plugged in.
    # The format is: Result = Numerator / (Term1 * Term2 * Term3 * Term4)
    print(f"{years_to_fade:.2f} years = {jnf_dose_bw1} / ({lux_level} * {hours_per_day} * {days_per_year} * {uv_damage_factor})")


calculate_fading_time()