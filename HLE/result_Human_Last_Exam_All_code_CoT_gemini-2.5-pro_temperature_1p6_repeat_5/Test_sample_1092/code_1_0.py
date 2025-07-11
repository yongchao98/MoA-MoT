def calculate_fade_time():
    """
    Calculates the time in years for a just noticeable fade to occur
    for a highly light-sensitive object.
    """

    # Step 1: Define the constants based on conservation science and the problem statement.
    # The total light exposure to cause a "Just Noticeable Fade" (JNF) for a
    # material with ISO Bluewool rating 1 is approximately 50,000 lux-hours.
    total_lux_hours_for_jnf = 50000

    # The light level on the object is 50 lux.
    illuminance_lux = 50

    # We assume a standard "daily" exposure of 8 hours for 365 days a year.
    hours_per_day = 8
    days_per_year = 365

    # Step 2: Calculate the total annual light exposure.
    annual_exposure = illuminance_lux * hours_per_day * days_per_year

    # Step 3: Calculate the time in years to reach the JNF threshold.
    years_to_fade = total_lux_hours_for_jnf / annual_exposure

    # Print the explanation and the final equation with all numbers.
    print(f"The calculation is based on the following formula:")
    print("Years to Fade = (Total Lux-Hours for JNF) / (Illuminance * Hours per Day * Days per Year)")
    print("\nPlugging in the values, the final equation is:")
    print(f"Years to Fade = {total_lux_hours_for_jnf} / ({illuminance_lux} * {hours_per_day} * {days_per_year})")
    print(f"\nThis results in approximately {years_to_fade:.2f} years for the next just noticeable fade to occur.")


if __name__ == '__main__':
    calculate_fade_time()