import math

def calculate_war_duration():
    """
    Calculates the duration of the Second Great War in StarCraft lore.

    The war started in 2504 and ended in 2506. This means the conflict
    spanned across three calendar years. To find the duration in years,
    rounded up, we count the number of years in the inclusive range.
    """
    start_year = 2504
    end_year = 2506

    # The duration spanning across calendar years is calculated as:
    # end_year - start_year + 1
    # This inherently handles the "rounded up" concept by counting any
    # part of a year as one full year in the span.
    duration = (end_year - start_year) + 1
    
    # We use math.ceil to demonstrate the "rounding up" from the raw difference,
    # although simply adding 1 to the difference of the years when counting
    # the number of involved years is more direct.
    # The war duration is slightly more than 2 years if we go by dates,
    # so ceiling(end - start) would not be accurate. The accurate way to
    # interpret "rounded up duration in years" when given start and end years
    # is to count the number of calendar years the event spanned.
    # So, 2504, 2505, 2506 -> 3 years.
    
    print(f"The Second Great War started in {start_year} and ended in {end_year}.")
    print(f"The number of calendar years spanned is {end_year} - {start_year} + 1 = {duration}")
    print(f"The duration of the war, rounded up to the nearest year, was {duration} years.")

calculate_war_duration()