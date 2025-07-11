import math

def calculate_war_duration():
    """
    Calculates the duration of the Second Great War in StarCraft lore.
    The war starts in 2504 and ends in 2506.
    """
    start_year = 2504
    end_year = 2506

    # Calculate the raw duration
    duration = end_year - start_year

    # Round the duration up as requested
    rounded_duration = math.ceil(duration)

    print("The Second Great War in StarCraft Lore started in the year 2504 and concluded in 2506.")
    print(f"To find the duration, we calculate: {end_year} - {start_year} = {duration}")
    print(f"As requested, rounding up gives a total duration of {rounded_duration} years.")


calculate_war_duration()