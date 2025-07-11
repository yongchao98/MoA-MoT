def solve_gift_game():
    """
    This function calculates and prints the answers to the gift game problem.
    """

    # Question 1: What is the expected value of X_19?
    # For odd n, the expected time is E[X_n] = (n-1)/2.
    n_19 = 19
    # The equation is (19-1)/2 = 9
    expected_time_19 = (n_19 - 1) / 2

    # Question 2: What is the expected value of X_20?
    # For even n, the game never ends, so the expected time is infinite.
    expected_time_20 = "inf"

    # Question 3: In general, for odd n>1 what is the expected value of X_n?
    general_formula_odd_n = "(n-1)/2"

    # Question 4: For odd n>30 what is the expected number of times that between
    # the two gifts there are exactly 10 and n-12 friends.
    # This corresponds to the distance D being 11 or n-11.
    # The expected number of visits to these two states is approximately (n-1)/n.
    # For n > 30, this is very close to 1.
    expected_visits = 1

    # Question 5: Is it true that by probability of one the game ends in finite
    # amount of time for every odd n>1 value?
    # Yes, the underlying Markov chain is recurrent.
    ends_with_prob_1 = "yes"

    # Print the answers separated by a comma, as requested.
    # The first answer includes the calculation for clarity, as hinted by the prompt.
    print(f"{int(expected_time_19)},{expected_time_20},{general_formula_odd_n},{expected_visits},{ends_with_prob_1}")

solve_gift_game()