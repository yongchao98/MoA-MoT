def solve_talbot_question():
    """
    This function calculates and prints the answers to the historical questions
    about the Talbot Settlement.
    """
    # Part 1: Number of settlers
    # Historical sources, such as the Dictionary of Canadian Biography, estimate the
    # population of the Talbot Settlement to be around 12,000 by 1823.
    # While the term "destitute" is subjective, many settlers were poor migrants.
    settlers_by_1823 = 12000

    # Part 2: Acreage difference
    # Colonel Talbot's original personal land grant in 1803 was 5,000 acres.
    original_grant_acres = 5000
    # He eventually oversaw the settlement of a much larger area,
    # estimated to be over 650,000 acres.
    total_administered_acres = 650000

    # Calculate the difference
    acreage_difference = total_administered_acres - original_grant_acres

    # Print the answers
    print(f"Between 1803 and 1823, approximately {settlers_by_1823} migrants settled in the Talbot Settlement.")
    print(f"The acreage Colonel Talbot eventually administered was {acreage_difference} acres larger than his original grant.")
    print(f"Calculation: {total_administered_acres} acres (total administered) - {original_grant_acres} acres (original grant) = {acreage_difference} acres.")

solve_talbot_question()