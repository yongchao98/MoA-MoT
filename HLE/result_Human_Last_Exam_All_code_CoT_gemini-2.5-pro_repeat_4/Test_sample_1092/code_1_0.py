import math

def calculate_fade_time():
    """
    Calculates the time in years for a Just Noticeable Fade (JNF) to occur
    for a highly light-sensitive material under specific lighting conditions.
    """
    # Step 1: Define constants based on the problem and conservation science standards.
    total_exposure_for_jnf_bw1 = 50000  # Lux-hours for a JNF on Bluewool 1 material.
    measured_illuminance_lux = 50       # Illuminance level in lux.
    uv_damage_factor = 3.0              # Assumption: UV-rich light is 3x more damaging.
    hours_per_day = 8.0                 # Assumption: 8 hours of exposure per day.
    days_per_year = 365.25              # Average days in a year to account for leap years.

    # Step 2: Calculate the effective illuminance.
    effective_illuminance_lux = measured_illuminance_lux * uv_damage_factor

    # Step 3: Calculate the total hours of exposure until JNF.
    total_hours_to_fade = total_exposure_for_jnf_bw1 / effective_illuminance_lux

    # Step 4: Convert total hours to total days.
    total_days_to_fade = total_hours_to_fade / hours_per_day

    # Step 5: Convert total days to years.
    total_years_to_fade = total_days_to_fade / days_per_year

    # --- Output the results step-by-step ---
    print("Calculation for Time to Just Noticeable Fade (JNF)")
    print("-" * 50)
    print(f"JNF Threshold for Bluewool 1: {int(total_exposure_for_jnf_bw1)} lux-hours")
    print(f"Measured Illuminance: {int(measured_illuminance_lux)} lux (UV-rich)")
    print("-" * 50)

    print("Equation: Years = (JNF Threshold / (Lux * UV Factor)) / Hours per Day / Days per Year\n")

    print(f"1. Effective Illuminance = {int(measured_illuminance_lux)} lux * {int(uv_damage_factor)} (UV factor) = {int(effective_illuminance_lux)} effective lux")
    print(f"2. Total Hours to Fade  = {int(total_exposure_for_jnf_bw1)} lux-hours / {int(effective_illuminance_lux)} effective lux = {total_hours_to_fade:.2f} hours")
    print(f"3. Total Days to Fade   = {total_hours_to_fade:.2f} hours / {int(hours_per_day)} hours/day = {total_days_to_fade:.2f} days")
    print(f"4. Total Years to Fade  = {total_days_to_fade:.2f} days / {days_per_year} days/year = {total_years_to_fade:.2f} years")
    print("-" * 50)
    print(f"Final Answer: The next just noticeable fade will occur in approximately {total_years_to_fade:.2f} years.")

if __name__ == '__main__':
    calculate_fade_time()
<<<0.11>>>