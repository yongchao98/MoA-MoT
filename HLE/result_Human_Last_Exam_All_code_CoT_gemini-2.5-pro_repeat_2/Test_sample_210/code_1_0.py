def solve_talbot_land_grant():
    """
    Calculates and prints answers to questions about the Talbot Settlement.
    """
    # Part 1: Number of migrants settled by 1823.
    # Historical sources state the population of the Talbot Settlement reached 12,000 by 1823.
    settlers_by_1823 = 12000

    # Part 2: Difference in acreage.
    # Colonel Talbot's original grant was 5,000 acres.
    original_grant_acres = 5000
    # He eventually administered or claimed over 650,000 acres across 28 townships.
    final_claimed_acres = 650000

    # Calculate how much larger the final acreage was.
    acreage_difference = final_claimed_acres - original_grant_acres

    # Print the answers.
    print(f"Number of migrants settled between 1803 and 1823: {settlers_by_1823}")
    print("\nCalculation for the difference in acreage:")
    print(f"The acreage he claimed ({final_claimed_acres}) was larger than the original grant ({original_grant_acres}).")
    print("The difference is:")
    print(f"{final_claimed_acres} - {original_grant_acres} = {acreage_difference}")
    print(f"\nThe acreage was {acreage_difference} acres larger.")

solve_talbot_land_grant()