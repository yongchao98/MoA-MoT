def solve_talbot_settlement_query():
    """
    This function calculates and prints answers to historical questions
    about the Talbot Settlement based on researched data.
    """
    
    # Part 1: Number of migrants by 1823
    # Historical sources estimate that by 1823, the population of the
    # Talbot Settlement had reached approximately 12,000 people.
    migrants_1823 = 12000

    # Part 2: Acreage difference
    # Colonel Talbot's original personal land grant in 1803 was 5,000 acres.
    original_grant_acres = 5000
    
    # Over the years, he came to administer the settlement of a vast territory,
    # estimated to be over 650,000 acres.
    total_claimed_acres = 650000
    
    # Calculate the difference
    acreage_difference = total_claimed_acres - original_grant_acres
    
    # Print the results
    print(f"Approximately {migrants_1823} migrants had settled in the Talbot Settlement between 1803 and 1823.")
    print("\nThe acreage he ultimately claimed was larger than his original grant by this amount:")
    print(f"{total_claimed_acres} (total acres) - {original_grant_acres} (original grant) = {acreage_difference} acres")

solve_talbot_settlement_query()