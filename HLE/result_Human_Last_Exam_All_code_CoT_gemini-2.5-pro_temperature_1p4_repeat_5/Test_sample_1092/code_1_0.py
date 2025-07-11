import math

def calculate_fade_time():
    """
    Calculates the time in years for a Just Noticeable Fade (JNF) to occur
    on a highly light-sensitive material based on conservation standards.
    """
    # --- Step 1: Define Known Values and Assumptions ---

    # Standard cumulative exposure in lux-hours to cause a JNF for a material
    # with ISO Bluewool rating of 1.
    lux_hours_to_jnf_standard = 50000

    # "UV-rich" light accelerates fading. A common assumption is that it
    # doubles the damage rate compared to standard filtered light.
    uv_damage_factor = 2

    # Light intensity the object is exposed to in lux.
    lux_level = 50

    # Assumed hours of exposure per day (based on a standard 8-hour museum day).
    hours_per_day = 8

    # Average number of days in a year to account for leap years.
    days_per_year = 365.25

    # --- Step 2: Adjust Fading Threshold for UV-Rich Light ---
    lux_hours_to_jnf_effective = lux_hours_to_jnf_standard / uv_damage_factor

    # --- Step 3: Calculate Total Annual Exposure ---
    annual_exposure = lux_level * hours_per_day * days_per_year

    # --- Step 4: Calculate Years to Fade ---
    years_to_jnf = lux_hours_to_jnf_effective / annual_exposure

    # --- Output the Explanation and Result ---
    print("This calculation determines the time until a Just Noticeable Fade (JNF) occurs.")
    print("-" * 60)
    print(f"JNF Threshold for Bluewool 1: {lux_hours_to_jnf_standard} lux-hours")
    print(f"UV-Rich Light Damage Factor: {uv_damage_factor}")
    print(f"Effective JNF Threshold: {int(lux_hours_to_jnf_effective)} lux-hours")
    print("-" * 60)

    print("The final calculation is:")
    print("Years to Fade = (Effective JNF Threshold) / (Annual Exposure)")
    print("Years to Fade = (JNF Threshold / UV Factor) / (Lux Level * Hours/Day * Days/Year)")
    print("\nPlugging in the numbers:")
    # The prompt requires printing each number in the final equation.
    print(f"Years = ({lux_hours_to_jnf_standard} / {uv_damage_factor}) / ({lux_level} * {hours_per_day} * {days_per_year})")
    print(f"Years = {int(lux_hours_to_jnf_effective)} / {int(annual_exposure)}")
    print(f"\nResult: It will take approximately {years_to_jnf:.3f} years for the next just noticeable fade to occur.")
    
    # Return the final numerical answer for the system.
    return years_to_jnf

# Execute the calculation and store the result.
final_answer = calculate_fade_time()
# The required output format is <<<answer>>>
# print(f"<<<{final_answer:.3f}>>>")