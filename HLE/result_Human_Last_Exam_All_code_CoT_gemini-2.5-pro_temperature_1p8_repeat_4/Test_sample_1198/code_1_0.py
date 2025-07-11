def solve_problem():
    """
    This function prints the answers to the two questions based on established mathematical results.
    """
    # For the first question: Is there any a>0 such that floor(a^n) = n mod 2 for every n>0?
    # The answer is Yes.
    answer_mod_2 = "Yes"

    # For the second question: Is there any a>0 such that floor(a^n) = n mod 3 for every n>0?
    # The answer is No.
    answer_mod_3 = "No"

    # Print the final answer in the required format.
    print(f"{answer_mod_2},{answer_mod_3}")

solve_problem()