import math

def calculate_war_duration():
    """
    Calculates the duration of the Second Great War in StarCraft lore.
    """
    # According to lore, the war starts in 2504 with Wings of Liberty.
    start_year = 2504

    # The war concludes in 2508 with the Into the Void epilogue campaign.
    end_year = 2508

    # Calculate the duration in years.
    duration = end_year - start_year

    # Round up the result as requested. Since the duration is a whole number,
    # rounding up does not change the value.
    final_duration = math.ceil(duration)

    # Print the equation with all numbers and the final answer.
    print(f"The Second Great War starts in the year {start_year} and ends in the year {end_year}.")
    print(f"Calculation: {end_year} - {start_year} = {duration}")
    print(f"The total duration of the war is {final_duration} years.")

if __name__ == "__main__":
    calculate_war_duration()