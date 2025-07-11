def solve_talbot_history():
    """
    Calculates and prints historical figures related to the Talbot Settlement.
    """
    # Historical data
    settlers_by_1823 = 12000
    original_grant_acres = 5000
    final_claimed_acres = 650000

    # Calculate the difference in acreage
    acreage_difference = final_claimed_acres - original_grant_acres

    # Print the answers
    print(f"Number of settlers between 1803 and 1823: {settlers_by_1823}")
    print("\nCalculation of the increase in acreage claimed:")
    print(f"{final_claimed_acres} (final acres) - {original_grant_acres} (original grant) = {acreage_difference} acres larger")

solve_talbot_history()