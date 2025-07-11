import math

def calculate_war_duration():
    """
    Calculates and prints the duration of the Second Great War in StarCraft lore.
    """
    # According to StarCraft lore, the Second Great War began in 2504
    # with the events of 'Wings of Liberty'.
    start_year = 2504

    # The war concluded in 2506 with the end of the 'Legacy of the Void' campaign.
    end_year = 2506

    # The war was active across the years 2504, 2505, and 2506.
    # While a direct subtraction (end_year - start_year) equals 2, the conflict
    # spanned parts of three calendar years, making its true duration greater than 2.
    # To find the duration in years rounded up, we can count the number of
    # calendar years the war was active in.
    duration_rounded_up = (end_year - start_year) + 1

    print(f"The Second Great War started in the year {start_year} and ended in {end_year}.")
    print("Since the conflict spanned across three calendar years, the duration rounded up is 3 years.")
    print("\nHere is the calculation:")
    print(f"{end_year} - {start_year} + 1 = {duration_rounded_up}")


calculate_war_duration()