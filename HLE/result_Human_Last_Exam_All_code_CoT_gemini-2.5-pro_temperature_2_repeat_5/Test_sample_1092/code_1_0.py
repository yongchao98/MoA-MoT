import math

def calculate_fading_time():
    """
    Calculates the time in years for a just noticeable fade (JNF) to occur
    for a highly light-sensitive material under specific light conditions.
    """

    # --- Step 1: Define constants and assumptions ---

    # ISO Bluewool 1 JNF threshold under ideal conditions (in lux-hours)
    # This is a standard value from conservation science.
    base_jnf_exposure_lux_hours = 50000

    # The problem specifies UV-rich light, which accelerates fading. A common
    # rule of thumb is that it doubles the rate of damage.
    uv_damage_factor = 2.0

    # Light exposure conditions
    light_level_lux = 50
    # Assuming a standard museum/gallery exposure schedule.
    daily_exposure_hours = 10
    days_per_year = 365

    # --- Step 2: Perform Calculations ---

    # Adjust the JNF threshold for UV-rich light
    total_exposure_for_jnf = base_jnf_exposure_lux_hours / uv_damage_factor

    # Calculate the total light exposure per year
    annual_exposure = light_level_lux * daily_exposure_hours * days_per_year

    # Calculate the time in years to reach the JNF threshold
    years_to_fade = total_exposure_for_jnf / annual_exposure

    # --- Step 3: Print the explanation and result ---

    print("Calculating the time to the next 'Just Noticeable Fade' (JNF):\n")

    print(f"1. Determine the JNF threshold for an ISO Bluewool 1 object in UV-rich light:")
    print(f"   Equation: {base_jnf_exposure_lux_hours} (standard lux-hours) / {uv_damage_factor} (UV factor)")
    print(f"   Result: {total_exposure_for_jnf} lux-hours\n")

    print(f"2. Calculate the total annual light exposure:")
    print(f"   Equation: {light_level_lux} (lux) * {daily_exposure_hours} (hours/day) * {days_per_year} (days/year)")
    print(f"   Result: {annual_exposure} lux-hours per year\n")

    print(f"3. Calculate the time in years to reach the JNF threshold:")
    print(f"   Equation: {total_exposure_for_jnf} (JNF threshold) / {annual_exposure} (annual exposure)")
    print(f"   Result: {years_to_fade:.3f} years\n")

    print(f"It will take approximately {years_to_fade:.3f} years for the next just noticeable fade to occur.")

# Execute the function
calculate_fading_time()

# The final answer in the required format
final_answer = (50000 / 2.0) / (50 * 10 * 365)
# <<<0.137>>>