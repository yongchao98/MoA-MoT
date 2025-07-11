def talbot_settlement_facts():
    """
    Calculates and prints facts about the Talbot Settlement between 1803 and 1823.
    """
    
    # Part 1: Number of settlers
    # Historical sources indicate that by 1823, around 20,000 people had settled in the area.
    settlers_by_1823 = 20000
    print(f"Between 1803 and 1823, an estimated {settlers_by_1823} migrants settled in the Talbot Settlement.")
    
    print("-" * 40)

    # Part 2: Acreage difference
    # Colonel Talbot's initial land grant was 5,000 acres.
    # He eventually came to administer a vast area of over 650,000 acres.
    original_grant_acres = 5000
    final_claimed_acres = 650000

    # Calculate the difference
    acreage_difference = final_claimed_acres - original_grant_acres

    print("To find how much larger the claimed acreage was than the original grant:")
    print(f"Final Acreage Claimed: {final_claimed_acres}")
    print(f"Original Land Grant: {original_grant_acres}")
    print("\nThe equation for the difference is:")
    print(f"{final_claimed_acres} (Final) - {original_grant_acres} (Original) = {acreage_difference} (Difference)")
    print(f"\nThe acreage he claimed was {acreage_difference} acres larger than the original land grant.")

talbot_settlement_facts()