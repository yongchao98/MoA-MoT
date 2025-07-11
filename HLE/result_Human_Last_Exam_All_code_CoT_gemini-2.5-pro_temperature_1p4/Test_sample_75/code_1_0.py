import math

def calculate_war_duration():
    """
    Calculates the duration of the Second Great War in StarCraft lore.
    The war starts in 2504 and ends in 2506.
    """
    start_year = 2504
    end_year = 2506

    # Calculate the raw duration
    raw_duration = end_year - start_year
    
    # Round up the result as per the request
    final_duration = math.ceil(raw_duration)

    print(f"The Second Great War started in the year {start_year} and ended in the year {end_year}.")
    print("To find the duration, we subtract the start year from the end year.")
    print(f"Equation: {end_year} - {start_year} = {raw_duration}")
    print(f"The total duration is {final_duration} years.")

calculate_war_duration()