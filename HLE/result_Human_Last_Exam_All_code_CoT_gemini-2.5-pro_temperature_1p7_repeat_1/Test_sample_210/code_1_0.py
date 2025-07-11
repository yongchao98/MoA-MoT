def solve_talbot_land_grant():
    """
    This function calculates and prints information about the Talbot Settlement.
    """
    # Part 1: Number of migrants
    # Historical sources estimate that by the early 1820s, the population
    # of the Talbot Settlement had grown to 20,000.
    migrant_settlers = 20000

    # Part 2: Acreage calculation
    # Colonel Talbot's original land grant in 1803.
    original_grant_acres = 5000

    # The total acreage he eventually claimed for his personal estate.
    final_claimed_acres = 65000

    # Calculate the difference.
    acreage_difference = final_claimed_acres - original_grant_acres

    # Print the full answer.
    print(f"An estimated {migrant_settlers} migrants settled in the Talbot Settlement between 1803 and 1823.")
    print(f"The acreage Colonel Talbot ultimately claimed ({final_claimed_acres} acres) was larger than the original land grant ({original_grant_acres} acres).")
    print(f"The difference is found by the equation: {final_claimed_acres} - {original_grant_acres} = {acreage_difference} acres.")

solve_talbot_land_grant()