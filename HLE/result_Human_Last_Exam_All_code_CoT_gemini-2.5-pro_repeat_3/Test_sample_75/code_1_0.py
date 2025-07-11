import math

def calculate_war_duration():
    """
    Calculates the duration of the Second Great War in StarCraft lore.
    The war starts with the Brood War (2500) and ends with the End War (2506).
    """
    start_year = 2500
    end_year = 2506

    # Calculate the raw duration
    duration = end_year - start_year

    # Round up the duration to the nearest whole year
    rounded_up_duration = math.ceil(duration)

    # Print the equation with the final rounded-up answer
    print(f"The Second Great War in StarCraft lore spans from {start_year} to {end_year}.")
    print("The duration in years, rounded up, is calculated as follows:")
    print(f"{end_year} - {start_year} = {rounded_up_duration}")

calculate_war_duration()