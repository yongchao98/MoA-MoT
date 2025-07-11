import math

def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur
    based on Bluewool rating, light level, and light type.
    """
    # Step 1: Define standard light dose for Bluewool 1
    lux_hours_for_bw1 = 50000  # Lux-hours for 1 JNF at BW1 with low-UV light

    # Step 2: Define adjustment factors
    # A damage factor of 2 is used for UV-rich light (e.g., daylight)
    uv_damage_factor = 2
    # Light intensity in lux
    light_intensity_lux = 50
    # Assuming a standard 8-hour "daily" exposure for an object on display
    hours_per_day = 8
    days_per_year = 365.25 # Using 365.25 to account for leap years

    # Step 3: Calculate the effective lux-hours needed for 1 JNF under UV-rich light
    effective_lux_hours_for_jnf = lux_hours_for_bw1 / uv_damage_factor

    # Step 4: Calculate the daily light dose
    daily_lux_hour_dose = light_intensity_lux * hours_per_day

    # Step 5: Calculate time to fade in years
    days_to_jnf = effective_lux_hours_for_jnf / daily_lux_hour_dose
    years_to_jnf = days_to_jnf / days_per_year

    # Output the final equation with all the numbers
    print("Equation to calculate years to next fade:")
    print(f"(({lux_hours_for_bw1} lux-hours / {uv_damage_factor} UV-factor) / ({light_intensity_lux} lux * {hours_per_day} hours/day)) / {days_per_year} days/year")
    print("\nCalculation:")
    print(f"({effective_lux_hours_for_jnf} / {daily_lux_hour_dose}) / {days_per_year} = {years_to_jnf:.3f} years")
    
    # Returning the final value for the formatted answer
    return years_to_jnf

# Run the calculation and store the result
final_answer = calculate_fade_time()

# The final answer is wrapped in <<<>>>
# Using .3f formatting to provide a reasonably precise answer.
print(f"\nIt will take approximately {final_answer:.3f} years for the next just noticeable fade to occur.")
print(f"<<<{final_answer:.3f}>>>")
