def solve_talbot_history():
    """
    This function calculates and prints answers to historical questions
    about the Talbot Settlement.
    """
    # Part 1: Number of migrants
    # Historical sources indicate that by the early 1820s, the population
    # of the Talbot Settlement had grown to approximately 20,000.
    migrants_settled = 20000

    # Part 2: Acreage difference
    # Colonel Talbot's original personal land grant was 5,000 acres.
    original_grant_acres = 5000
    # Through his role as a land settlement commissioner, he eventually
    # controlled the settlement of about 650,000 acres.
    total_claimed_acres = 650000

    # Calculate the difference
    acreage_difference = total_claimed_acres - original_grant_acres

    # Print the results
    print(f"Approximately {migrants_settled} destitute migrants settled as a result of the Talbot land grant between 1803 and 1823.")
    print("\nThe acreage Colonel Talbot eventually claimed was larger than his original personal grant by:")
    print(f"{total_claimed_acres} - {original_grant_acres} = {acreage_difference} acres.")

solve_talbot_history()