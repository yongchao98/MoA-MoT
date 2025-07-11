def solve_talbot_land_grant():
    """
    Calculates and prints information about the Talbot Settlement.
    This includes the number of settlers by 1823 and the difference
    in land acreage from the original grant to the final claim.
    """
    # Historical data points
    settlers_by_1823 = 12000
    original_grant_acres = 5000
    final_claimed_acres = 650000

    # Calculate the difference in acreage
    acreage_difference = final_claimed_acres - original_grant_acres

    # Print the answer for the number of settlers
    print(f"It is estimated that {settlers_by_1823} migrants had settled in the Talbot Settlement between 1803 and 1823.")
    print("\nTo find how many acres larger the final claimed acreage was than the original grant, we perform the following subtraction:")
    
    # Print the equation with all numbers, as requested
    print(f"Final Claimed Acreage ({final_claimed_acres}) - Original Grant Acreage ({original_grant_acres}) = Difference ({acreage_difference})")
    
    # Print the final answer for the acreage difference
    print(f"\nThe acreage Colonel Talbot eventually claimed was {acreage_difference} acres larger than his original land grant.")

# Execute the function
solve_talbot_land_grant()