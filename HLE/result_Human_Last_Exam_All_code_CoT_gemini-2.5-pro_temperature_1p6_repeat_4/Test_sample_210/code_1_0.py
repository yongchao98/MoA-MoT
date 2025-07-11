def solve_talbot_land_grant():
    """
    Calculates and prints the answers to historical questions about the Talbot Settlement.
    """
    # Part 1: Number of settlers
    # Historical sources estimate that by the early 1820s, the population of the
    # Talbot Settlement had grown to approximately 20,000 people.
    settlers_by_1823 = 20000

    # Part 2: Acreage calculation
    # Colonel Talbot's initial land grant was 5,000 acres.
    initial_grant_acres = 5000
    # Through his role as a settlement agent, he amassed a personal estate of approximately 65,000 acres.
    final_claimed_acres = 65000

    # Calculate the difference
    acreage_difference = final_claimed_acres - initial_grant_acres

    # Print the results
    print(f"By 1823, approximately {settlers_by_1823} destitute migrants had settled in the Talbot Settlement.")
    print("To find how much larger the claimed acreage was than the original grant, we calculate the difference:")
    print(f"Final Claimed Acres - Initial Grant Acres = Difference")
    print(f"{final_claimed_acres} acres - {initial_grant_acres} acres = {acreage_difference} acres")
    print(f"\nThe acreage he claimed was {acreage_difference} acres larger than the original grant.")

solve_talbot_land_grant()

# <<<20,000 migrants settled, and the claimed acreage was 60,000 acres larger.>>>