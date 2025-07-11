def solve_talbot_settlement():
    """
    This function calculates and prints information about the Talbot Settlement.
    """
    # Historical data points for the Talbot Settlement.
    # The original land grant to Colonel Talbot was 5,000 acres.
    original_grant_acres = 5000
    
    # By the end of the 1820s, the population of the settlement he oversaw was around 20,000.
    settlers_by_1823 = 20000
    
    # Colonel Talbot eventually became responsible for settling over half a million acres,
    # with a commonly cited figure being around 650,000 acres across 29 townships.
    total_managed_acres = 650000

    # Calculate the difference in acreage
    acreage_difference = total_managed_acres - original_grant_acres

    # Print the answer for the number of migrants
    print(f"Number of destitute migrants who settled between 1803 and 1823: {settlers_by_1823}")

    # Print the answer for the acreage difference, showing the equation
    print("\nThe acreage he ultimately claimed was larger than the original land grant by:")
    print(f"{total_managed_acres} acres - {original_grant_acres} acres = {acreage_difference} acres")

solve_talbot_settlement()