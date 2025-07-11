def solve_talbot_settlement_question():
    """
    This function provides information about the Talbot Settlement and calculates
    the difference in land acreage claimed by Colonel Talbot.
    """
    
    # Part 1: Number of settlers
    # Historical sources indicate that by the early 1820s, the population
    # of the Talbot Settlement had grown to approximately 20,000 people.
    settlers_by_1823 = 20000

    # Part 2: Acreage calculation
    # Colonel Talbot's original personal land grant in 1803 was 5,000 acres.
    original_grant_acreage = 5000
    
    # Over his career, through his influential position and further grants,
    # he personally acquired a total of about 60,000 acres.
    total_claimed_acreage = 60000

    # Calculate the difference
    acreage_difference = total_claimed_acreage - original_grant_acreage

    # Print the answers
    print(f"As a result of the 1803 land grant, approximately {settlers_by_1823:,} destitute migrants had settled in the Talbot Settlement by 1823.")
    
    print("\nThe acreage Colonel Talbot eventually claimed was significantly larger than his original grant.")
    print("To find the difference, we perform the following subtraction:")
    print(f"{total_claimed_acreage} (total acres claimed) - {original_grant_acreage} (original grant acres) = {acreage_difference} (difference in acres)")
    
    print(f"\nThis means the acreage he claimed was {acreage_difference:,} acres larger than the original land grant.")

solve_talbot_settlement_question()