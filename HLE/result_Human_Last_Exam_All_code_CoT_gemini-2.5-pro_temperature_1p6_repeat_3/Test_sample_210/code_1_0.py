def solve_talbot_land_grant():
    """
    Calculates and prints information about the Talbot Settlement.
    """
    # Historical data
    # Sources indicate the population of the Talbot Settlement was about 20,000 by 1823.
    settlers_by_1823 = 20000
    
    # Colonel Talbot's original land grant was 5,000 acres.
    original_grant_acres = 5000
    
    # He eventually administered the settlement of over 650,000 acres.
    final_claimed_acres = 650000

    # Calculate the difference in acreage
    acreage_difference = final_claimed_acres - original_grant_acres

    # Print the answer to the first question
    print(f"Between 1803 and 1823, it is estimated that {settlers_by_1823} migrants settled in the Talbot Settlement.")
    print("-" * 20)

    # Print the answer to the second question, showing the calculation
    print("To find how much larger the claimed acreage was than the original grant:")
    print(f"Final Claimed Acreage ({final_claimed_acres}) - Original Grant ({original_grant_acres}) = Acreage Difference ({acreage_difference})")
    print(f"\nThe acreage Colonel Talbot claimed was {acreage_difference:,} acres larger than his original land grant.")

solve_talbot_land_grant()