def solve_talbot_settlement_questions():
    """
    Calculates and prints answers to questions about the Talbot Settlement.
    """
    # Historical data for the period 1803-1823
    migrants_settled = 20000
    original_land_grant_acres = 5000
    total_claimed_acres = 650000

    # Calculate the difference in acreage
    acreage_difference = total_claimed_acres - original_land_grant_acres

    # Print the answer for the number of migrants
    print(f"Approximately {migrants_settled} destitute migrants settled in the Talbot Settlement between 1803 and 1823.")
    
    # Print the answer for the acreage difference, showing the full equation
    print("\nThe acreage Colonel Talbot eventually claimed was larger than his original grant.")
    print(f"The calculation is: {total_claimed_acres} (total acres) - {original_land_grant_acres} (original grant) = {acreage_difference} acres.")
    print(f"This means the acreage he claimed was {acreage_difference} acres larger.")

# Execute the function
solve_talbot_settlement_questions()