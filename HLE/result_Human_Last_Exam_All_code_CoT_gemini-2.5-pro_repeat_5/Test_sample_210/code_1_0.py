def solve_talbot_settlement_query():
    """
    Calculates and prints information about the Talbot Settlement.
    """
    # Historical data based on research
    migrants_by_1823 = 20000
    original_grant_acres = 5000
    final_administered_acres = 650000

    # Calculate the difference in acreage
    acreage_difference = final_administered_acres - original_grant_acres

    # Print the answer to the first question
    print(f"Approximately {migrants_by_1823} migrants settled in the Talbot Settlement between 1803 and 1823.")
    print("-" * 20)

    # Print the answer to the second question, showing the equation
    print("The acreage Colonel Talbot eventually administered was larger than his original grant.")
    print("The calculation is as follows:")
    print(f"Final Acreage ({final_administered_acres}) - Original Grant ({original_grant_acres}) = Difference ({acreage_difference})")
    print(f"\nThis means the land he administered was {acreage_difference} acres larger than his original grant.")

solve_talbot_settlement_query()