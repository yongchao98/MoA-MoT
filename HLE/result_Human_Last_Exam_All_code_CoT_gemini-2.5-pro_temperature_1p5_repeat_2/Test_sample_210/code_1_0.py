def solve_talbot_history():
    """
    This function calculates and prints answers to historical questions
    about Colonel Thomas Talbot's land grant and settlement.
    """
    # Part 1: Number of settlers
    # Historical data shows that by 1823, the population of the Talbot
    # Settlement was approximately 12,000.
    start_year = 1803
    end_year = 1823
    settler_population = 12000

    print(f"Between {start_year} and {end_year}, the number of destitute migrants who settled was approximately {settler_population:,}.")
    print("-" * 40)

    # Part 2: Acreage difference
    # Historical data on the original grant and the eventually claimed/administered acreage.
    original_grant_acres = 5000
    claimed_acres = 650000

    # Calculate the difference
    acreage_difference = claimed_acres - original_grant_acres

    print("To find how much larger the claimed acreage was than the original grant:")
    print(f"Equation: {claimed_acres:,} (claimed acres) - {original_grant_acres:,} (original grant acres) = {acreage_difference:,} acres")
    print(f"The acreage he eventually claimed was {acreage_difference:,} acres larger than the original grant.")

solve_talbot_history()