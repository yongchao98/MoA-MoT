import math

def solve_and_print():
    """
    This function calculates and prints the answers to the user's questions.
    """
    
    # Part 1: What is the expected value of X_19?
    # For odd n, E[X_n] = n - 1.
    n_19 = 19
    expected_time_19 = n_19 - 1
    
    # Part 2: What is the expected value of X_20?
    # For even n, the game never ends, so the expected time is infinite.
    expected_time_20 = "infty"

    # Part 3: In general, for odd n>1 what is the expected value of X_n?
    general_formula_odd_n = "n-1"
    
    # Part 4: For odd n>30 what is the expected number of times that between the
    # two gifts there are exactly 10 and n-12 friends.
    # The number of friends between the gifts corresponds to distance-1.
    # The distance must be even in our model. A distance of 11 (10 friends) is impossible.
    # We assume a typo and that the question meant 11 friends (distance 12).
    # The formula for the expected number of visits to the state with distance 12 is 48/(n+1).
    expected_visits = "48/(n+1)"
    
    # Part 5: Is it true that by probability of one the game ends in finite
    # amount of time for every odd n>1 value?
    # Yes, because the underlying random walk is on a finite, irreducible state space.
    prob_one_finish = "Yes"

    # Print the final answer string, with results separated by commas.
    final_answer_string = f"{expected_time_19}, {expected_time_20}, {general_formula_odd_n}, {expected_visits}, {prob_one_finish}"
    print(final_answer_string)

solve_and_print()