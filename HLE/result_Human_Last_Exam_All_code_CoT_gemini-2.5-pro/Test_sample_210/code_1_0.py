def solve_talbot_land_grant():
    """
    This function calculates and prints information about the
    Talbot Settlement between 1803 and 1823.
    """
    # Historical data based on research.
    # By 1823, Talbot had overseen the settlement of 20,000 people.
    number_of_migrants = 20000
    start_year = 1803
    end_year = 1823

    # Colonel Talbot's original land grant was 5,000 acres.
    original_grant_acres = 5000

    # He eventually acquired personal title to approximately 65,000 acres.
    final_claimed_acres = 65000

    # Calculate how many more acres he claimed than his original grant.
    acreage_difference = final_claimed_acres - original_grant_acres

    # Print the answer to the first question.
    print(f"As a result of the land grant, {number_of_migrants} destitute migrants settled between {start_year} and {end_year}.")

    # Print the answer to the second question, showing the equation.
    print("\nThe acreage he ultimately claimed was larger than the original grant by this many acres:")
    print(f"{final_claimed_acres} (Final Claimed) - {original_grant_acres} (Original Grant) = {acreage_difference} acres")

solve_talbot_land_grant()