def solve():
    """
    This function calculates and prints the answers to the user's questions.
    """
    # Question 1: What is the expected value of X_19?
    # For an odd n, the expected time is given by the formula (n^2 - 1) / 2.
    n_19 = 19
    x_19_numerator = n_19**2 - 1
    x_19_denominator = 2
    x_19_val = x_19_numerator / x_19_denominator

    # Question 2: And for X_20?
    # For an even n, the distance between gifts is always odd, but the end-condition
    # requires an even distance (2 or n-2). Thus, the game never ends.
    x_20_val = "inf"

    # Question 3: In general, for odd n>1 what is the expected value of X_n?
    # This is the general formula derived above.
    x_n_formula = "(n**2 - 1) / 2"

    # Question 4: For odd n>30 what is the expected number of times that between
    # the two gifts there are exactly 10 and n-12 friends?
    # This corresponds to a distance of 10+1=11. The expected number of visits to
    # the state with distance 11, starting from distance 1, is n-11.
    visits_formula = "n - 11"

    # Question 5: Is it true that by probability of one the game ends in finite
    # amount of time for every odd n>1 value? Answer with yes or no.
    # Yes, for odd n, the game is an absorbing Markov chain where absorption
    # is certain.
    prob_one_val = "yes"

    # Print all answers in a comma-separated format.
    print(f"{int(x_19_val)},{x_20_val},{x_n_formula},{visits_formula},{prob_one_val}")

solve()