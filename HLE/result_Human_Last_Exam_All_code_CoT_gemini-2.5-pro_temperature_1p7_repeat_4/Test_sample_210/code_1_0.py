def solve_talbot_questions():
    """
    This function provides answers to historical questions about the Talbot Settlement.
    It prints the number of settlers by 1823 and calculates the difference
    between Talbot's total controlled land and his original grant.
    """
    
    # --- Historical Data ---
    # The approximate number of settlers in the Talbot Settlement by 1823.
    # While the term "destitute migrants" is specific, historical sources typically
    # refer to the total population, which reached this number.
    settlers_by_1823 = 20000

    # The size of Colonel Talbot's original land grant in 1803 in acres.
    original_land_grant_acres = 5000

    # The total acreage Talbot eventually controlled or administered the settlement of.
    total_controlled_acres = 650000

    # --- Calculation ---
    # Calculate how many more acres he controlled compared to his initial grant.
    acreage_difference = total_controlled_acres - original_land_grant_acres

    # --- Output ---
    # Answer the first question regarding the number of settlers.
    print(f"1. Approximately {settlers_by_1823} migrants had settled in the Talbot Settlement between 1803 and 1823.")
    print("")

    # Answer the second question regarding the acreage difference, showing the equation.
    print("2. The acreage he claimed was larger than the original land grant by the following amount:")
    print(f"{total_controlled_acres} acres (Total Claimed) - {original_land_grant_acres} acres (Original Grant) = {acreage_difference} acres")

solve_talbot_questions()