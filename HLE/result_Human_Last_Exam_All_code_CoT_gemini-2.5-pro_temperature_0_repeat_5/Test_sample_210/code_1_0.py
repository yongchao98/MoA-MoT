def solve_talbot_history():
    """
    This function provides historical figures related to the Talbot Settlement
    and calculates the difference in land acreage.
    """
    # By 1823, the population of the Talbot Settlement is estimated to have been around 12,000.
    settlers_by_1823 = 12000

    # Colonel Talbot's original personal land grant in 1803 was 5,000 acres.
    original_grant_acres = 5000

    # He eventually came to administer and control a vast area of over 650,000 acres.
    total_claimed_acres = 650000

    # Calculate the difference in acreage.
    acreage_difference = total_claimed_acres - original_grant_acres

    print(f"An estimated {settlers_by_1823} destitute migrants settled in the Talbot Settlement between 1803 and 1823.")
    print("\nTo find the difference in acreage:")
    print(f"The acreage he eventually claimed was {total_claimed_acres} acres.")
    print(f"The original land grant was {original_grant_acres} acres.")
    print(f"The claimed acreage was {total_claimed_acres} - {original_grant_acres} = {acreage_difference} acres larger than the original grant.")

solve_talbot_history()