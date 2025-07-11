def solve_talbot_settlement_query():
    """
    Calculates and prints answers to questions about the Talbot Settlement.
    """

    # Data based on historical records
    migrants_by_1823 = 12000
    original_grant_acres = 5000
    total_claimed_acres = 650000

    # Calculate the difference in acreage
    acreage_difference = total_claimed_acres - original_grant_acres

    # --- Print the answers ---

    # Answer to the first question
    print(f"Number of destitute migrants settled by 1823: {migrants_by_1823}")

    # Answer to the second question
    print("\nCalculation for the difference in acreage:")
    print(f"The acreage Colonel Talbot eventually claimed was {total_claimed_acres:,} acres.")
    print(f"The original land grant was {original_grant_acres:,} acres.")
    print(f"The difference is: {total_claimed_acres:,} - {original_grant_acres:,} = {acreage_difference:,} acres.")
    print(f"\nThe claimed acreage was {acreage_difference:,} acres larger than the original grant.")

solve_talbot_settlement_query()