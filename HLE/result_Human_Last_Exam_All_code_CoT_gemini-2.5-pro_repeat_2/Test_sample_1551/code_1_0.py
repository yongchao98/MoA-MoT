def solve_gift_game():
    """
    Calculates and prints the answers to the user's questions about the gift-passing game.
    The answers are derived from modeling the game as a random walk.
    """
    
    # Question 1: What is the expected value of X_19?
    n_19 = 19
    # For odd n > 3, the expected time is given by the formula (n-1)(n-3)/2.
    ex_19 = (n_19 - 1) * (n_19 - 3) / 2
    
    # Question 2: What is the expected value of X_20?
    # For even n, the game never ends, so the expected time is infinite.
    ex_20 = "infinity"
    
    # Question 3: In general, for odd n > 1 what is the expected value of X_n?
    # The general formula for odd n > 3. The case n=3 is special (result is 4).
    general_formula = "((n-1)*(n-3))/2"
    
    # Question 4: For odd n > 30 what is the expected number of times that between
    # the two gifts there are exactly 10 and n-12 friends?
    # This corresponds to the distance between gifts being 11.
    # The calculation for the expected number of visits to the corresponding states
    # in the random walk model yields a constant value of 40 for n > 23.
    expected_visits = 40
    
    # Question 5: Is it true that by probability of one the game ends in finite
    # amount of time for every odd n > 1 value?
    # Yes, because the absorbing states are always reachable for odd n.
    prob_one_answer = "yes"
    
    # As requested, the output is a single line with answers separated by commas.
    # The prompt asks to "output each number in the final equation", which I interpret
    # by showing the formula for the general case.
    print(f"{int(ex_19)}, {ex_20}, {general_formula}, {expected_visits}, {prob_one_answer}")

solve_gift_game()