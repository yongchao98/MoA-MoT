def solve_talbot_settlement_query():
    """
    This function provides information about the Talbot Settlement
    between 1803 and 1823 based on historical data.
    """

    # Historical data for the Talbot Settlement
    # By the early 1820s (specifically 1822-1823), the population of the settlement had grown significantly.
    settlers_by_1823 = 20000

    # The original land grant to Colonel Talbot in 1803.
    original_grant_acres = 5000

    # Over the years, Talbot became the administrator for a much larger territory.
    total_claimed_acres = 668000

    # Calculate the difference in acreage
    acreage_difference = total_claimed_acres - original_grant_acres

    # Print the answer to the first question
    print(f"Between 1803 and 1823, approximately {settlers_by_1823} migrants settled in the Talbot Settlement.")

    # Print the answer to the second question, including the equation
    print("\nThe acreage Colonel Talbot eventually claimed was significantly larger than his original grant.")
    print(f"The calculation for the difference is:")
    print(f"{total_claimed_acres} (total acres claimed) - {original_grant_acres} (original grant acres) = {acreage_difference} acres")

solve_talbot_settlement_query()