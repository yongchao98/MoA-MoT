def solve_talbot_questions():
    """
    This function provides answers to questions about Colonel Thomas Talbot's land grant
    based on historical data.
    """
    
    # Number of destitute migrants who settled between 1803 and 1823.
    # By 1823, historical records indicate he had placed 20,000 people on the land.
    migrants_settled = 20000

    # The original land grant Talbot received in 1803 was 5,000 acres.
    original_grant_acres = 5000

    # The total acreage he eventually came to administer was over 650,000 acres.
    final_claimed_acres = 650000

    # Calculate how much larger the final claimed acreage was than the original grant.
    acreage_difference = final_claimed_acres - original_grant_acres

    # Print the answer regarding the number of settlers.
    print(f"Between 1803 and 1823, {migrants_settled} destitute migrants settled in the Talbot Settlement.")
    
    # Print the answer regarding the acreage difference, showing the full equation.
    print("\nThe acreage he ultimately claimed was larger than the original grant by the following amount:")
    print(f"{final_claimed_acres} - {original_grant_acres} = {acreage_difference} acres")

solve_talbot_questions()