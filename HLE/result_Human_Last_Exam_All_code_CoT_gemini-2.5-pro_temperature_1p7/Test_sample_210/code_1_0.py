def solve_talbot_settlement_query():
    """
    Calculates and prints information about the Talbot Settlement.
    """
    # Historical data based on research.
    num_migrants_by_1823 = 20000
    original_grant_acres = 5000
    claimed_acres = 650000

    # Calculate the difference in acreage.
    acreage_difference = claimed_acres - original_grant_acres

    # Print the answer to the first question.
    print(f"Between 1803 and 1823, approximately {num_migrants_by_1823} migrants settled as a result of the Talbot land grant.")

    # Print the answer and the calculation for the second question.
    print("\nTo find how much larger the claimed acreage was than the original grant:")
    print(f"Claimed Acreage: {claimed_acres}")
    print(f"Original Grant: {original_grant_acres}")
    print(f"The calculation is: {claimed_acres} - {original_grant_acres} = {acreage_difference}")
    print(f"The acreage he claimed was {acreage_difference} acres larger than the original land grant.")

# Execute the function to display the results.
solve_talbot_settlement_query()