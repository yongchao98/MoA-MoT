def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur
    on a highly light-sensitive material.
    """
    # Step 1: Define the constants based on conservation science standards.
    # Total light exposure (in lux-hours) to cause a Just Noticeable Fade (JNF)
    # for a material with ISO Bluewool rating 1.
    total_jnf_exposure = 50000

    # Light intensity in lux.
    light_intensity_lux = 50

    # Assuming a standard museum exhibition schedule for "daily" exposure.
    hours_per_day = 8
    days_per_year = 365

    # Step 2: Calculate the annual light exposure rate.
    annual_exposure = light_intensity_lux * hours_per_day * days_per_year

    # Step 3: Calculate the time in years to reach the JNF threshold.
    years_to_fade = total_jnf_exposure / annual_exposure

    # Print the explanation and the final equation with the result.
    print(f"To find the time until a just noticeable fade, we divide the total required exposure by the annual exposure rate.")
    print(f"The total exposure threshold for a Bluewool 1 material is {total_jnf_exposure} lux-hours.")
    print(f"The annual exposure is calculated from {light_intensity_lux} lux, for {hours_per_day} hours a day, {days_per_year} days a year.\n")
    print("The final equation is:")
    print(f"{total_jnf_exposure} / ({light_intensity_lux} * {hours_per_day} * {days_per_year}) = {years_to_fade:.3f} years")

if __name__ == "__main__":
    calculate_fade_time()
    # The final numerical answer is approximately 0.342
    # years_to_fade = 50000 / (50 * 8 * 365) = 50000 / 146000 = 0.34246...
    final_answer = 50000 / (50 * 8 * 365)
    # print(f"\n<<<{final_answer:.3f}>>>")