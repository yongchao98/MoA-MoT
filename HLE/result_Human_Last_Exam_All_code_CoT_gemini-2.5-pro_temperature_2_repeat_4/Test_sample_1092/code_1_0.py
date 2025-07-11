import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur on a
    highly light-sensitive object under specific lighting conditions.
    """
    
    # --- Constants and Assumptions ---

    # JNF (Just Noticeable Fade) threshold for ISO Bluewool 1 materials in lux-hours.
    # This is a standard value in conservation science for highly fugitive materials.
    jnf_threshold_lux_hours = 50000
    
    # Given light level in lux.
    illuminance_lux = 50
    
    # Damage factor for "UV-rich" light. Standard incandescent light has a factor of 1.
    # UV-rich sources like unfiltered daylight or fluorescent lights are more damaging. A factor of 3 is a reasonable estimate.
    uv_damage_factor = 3
    
    # Assumed daily exposure time in a museum setting.
    hours_per_day = 8
    
    # Number of days in a year.
    days_per_year = 365

    # --- Calculation ---
    
    # 1. Calculate the total annual effective light exposure.
    # This is the damaging equivalent of the light dose the object receives in one year.
    annual_effective_exposure = illuminance_lux * uv_damage_factor * hours_per_day * days_per_year
    
    # 2. Calculate the number of years to reach the JNF threshold.
    years_to_jnf = jnf_threshold_lux_hours / annual_effective_exposure
    
    # --- Output ---
    
    print(f"Based on the given conditions and standard conservation assumptions:")
    print(f"An object with ISO Bluewool Rating 1 will experience a just noticeable fade in approximately {years_to_jnf:.2f} years.")
    print("\nCalculation breakdown:")
    print(f"Years = {jnf_threshold_lux_hours} / ({illuminance_lux} lux * {uv_damage_factor} [UV factor] * {hours_per_day} hours/day * {days_per_year} days/year)")

calculate_fade_time()
<<<0.11>>>