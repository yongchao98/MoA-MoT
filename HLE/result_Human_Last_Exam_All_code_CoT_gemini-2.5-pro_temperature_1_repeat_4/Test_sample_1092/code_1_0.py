import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur on a
    highly sensitive material under specific light conditions.
    """
    # --- Step 1: Define constants and inputs from the problem ---
    light_intensity = 50  # in lux
    iso_bluewool_rating = 1
    
    # Standard total exposure (in lux-hours) for a Just Noticeable Fade (JNF)
    # for ISO Bluewool 1. This is a standard figure in conservation science.
    jnf_total_exposure = 50000

    # --- Step 2: Define assumptions based on the problem description ---
    # "UV-rich light" is more damaging. A common assumption is a 3x damage factor
    # compared to standard low-UV gallery lighting.
    uv_damage_factor = 3
    
    # "Daily" exposure is assumed to be 8 hours, a typical museum/gallery day.
    hours_per_day = 8
    
    # Using 365.25 to account for leap years on average.
    days_per_year = 365.25
    
    # --- Step 3 & 4: Calculate the total effective light exposure dose per year ---
    # This is the total amount of damaging light the object receives in one year.
    # Formula: lux * hours/day * days/year * damage_factor
    effective_annual_exposure = light_intensity * hours_per_day * days_per_year * uv_damage_factor

    # --- Step 5: Calculate the number of years to reach the JNF ---
    # Formula: Total exposure for JNF / Effective annual exposure
    years_to_fade = jnf_total_exposure / effective_annual_exposure

    # --- Output the results ---
    print("This calculation determines the time to a 'just noticeable fade' based on conservation science principles.")
    print("The formula is: Years = Total JNF Exposure / (Light Intensity * Hours/Day * Days/Year * UV Factor)\n")
    print("Here are the values used in the calculation:")
    print(f"- Total JNF Exposure for ISO Bluewool {iso_bluewool_rating}: {jnf_total_exposure:,} lux-hours")
    print(f"- Light Intensity: {light_intensity} lux")
    print(f"- Daily Exposure Time (Assumption): {hours_per_day} hours/day")
    print(f"- UV Damage Factor for 'UV-rich' light (Assumption): {uv_damage_factor}\n")

    print("Final Equation:")
    # We print each number as requested by the user.
    print(f"{jnf_total_exposure} / ({light_intensity} * {hours_per_day} * {days_per_year} * {uv_damage_factor}) = {years_to_fade:.3f} years")

calculate_fade_time()
<<<0.114>>>