import math

def calculate_fade_time():
    """
    Calculates the time in years for a highly light-sensitive object
    to experience a "Just Noticeable Fade" (JNF) under specific lighting.
    """

    # --- Step 1: Define Constants and Assumptions ---
    # According to conservation science, the exposure to cause one JNF for a
    # material with ISO Bluewool Rating 1 is approximately 20,000 lux-hours.
    # This value accounts for a standard UV-rich light source.
    jnf_threshold_lux_hours = 20000

    # The light intensity the object is exposed to.
    light_intensity_lux = 50

    # We assume "daily" exposure in a museum context is 8 hours per day.
    hours_per_day = 8

    # We calculate the total exposure over a full year.
    days_per_year = 365

    print("This script calculates the time until the next 'Just Noticeable Fade' (JNF) occurs.")
    print(f"The JNF threshold for a highly sensitive material (ISO Bluewool 1) is {jnf_threshold_lux_hours} lux-hours.")
    print("-" * 70)

    # --- Step 2: Calculate Annual Light Exposure ---
    print("Step 1: Calculate the total light exposure per year.")
    annual_exposure = light_intensity_lux * hours_per_day * days_per_year
    print(f"Annual Exposure = {light_intensity_lux} lux * {hours_per_day} hours/day * {days_per_year} days/year = {annual_exposure} lux-hours per year.")
    print("-" * 70)

    # --- Step 3: Calculate Years to Fade ---
    print("Step 2: Calculate the number of years to reach the JNF threshold.")
    years_to_fade = jnf_threshold_lux_hours / annual_exposure
    
    # As requested, here is the final equation showing each number:
    print("The final equation is:")
    print(f"Years to Fade = {jnf_threshold_lux_hours} / ({light_intensity_lux} * {hours_per_day} * {days_per_year})")
    
    print("-" * 70)
    print(f"Result: It will take approximately {years_to_fade:.2f} years for the next just noticeable fade to occur.")
    
    # Return the final numerical answer for the submission format
    return round(years_to_fade, 2)

# Execute the function and get the final answer
final_answer = calculate_fade_time()

# The final answer is enclosed in <<< >>>
print(f"\n<<< {final_answer} >>>")