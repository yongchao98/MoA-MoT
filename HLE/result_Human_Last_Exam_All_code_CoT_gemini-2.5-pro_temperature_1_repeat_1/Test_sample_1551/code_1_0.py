def solve_and_print():
    """
    This function calculates and prints the answers to the user's questions based on the derivations above.
    """

    # --- Question 1: Expected value of X_19 ---
    n_for_X19 = 19
    # The general formula for odd n is (n-1)/2.
    # For n=19, the equation is (19 - 1) / 2.
    # The numbers in the equation are 19, 1, and 2.
    numerator = n_for_X19 - 1
    denominator = 2
    e_x19_val = int(numerator / denominator)
    ans_1 = str(e_x19_val)

    # --- Question 2: Expected value of X_20 ---
    # For even n, the game never ends, so the expected time is infinite.
    ans_2 = "infinity"

    # --- Question 3: General formula for E[X_n] for odd n > 1 ---
    # The formula is (n-1)/2. The numbers in the equation are 1 and 2.
    num_in_formula_1 = 1
    num_in_formula_2 = 2
    ans_3 = f"(n-{num_in_formula_1})/{num_in_formula_2}"

    # --- Question 4: Expected number of visits for odd n > 30 ---
    # The formula is 2*(n-11)/n. The numbers in the equation are 2 and 11.
    num_in_formula_2_b = 2
    num_in_formula_11 = 11
    ans_4 = f"{num_in_formula_2_b}*(n-{num_in_formula_11})/n"

    # --- Question 5: Does the game end with probability 1 for odd n? ---
    # Yes, the random walk is guaranteed to be absorbed.
    ans_5 = "Yes"

    # --- Final Output ---
    # Combine all answers into a single comma-separated string as requested.
    final_output = f"{ans_1}, {ans_2}, {ans_3}, {ans_4}, {ans_5}"
    print(final_output)

solve_and_print()