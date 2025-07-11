import math

def solve_and_print():
    """
    This function calculates the answers to the user's questions and prints them
    in a comma-separated format.
    """

    # Question 1: What is the expected value of X_19?
    # For odd n, the formula is E[X_n] = (n^2 - 1) / 2.
    n_19 = 19
    # The numbers in the equation are n_19, 2, and 1.
    # The calculation is (19^2 - 1) / 2.
    ex_19 = (n_19**2 - 1) // 2

    # Question 2: What is the expected value of X_20?
    # For even n, the difference in gift positions starts odd and can only change by
    # an even number, so it can never become 0. The game never ends.
    ex_20 = "∞"

    # Question 3: In general, for odd n>1 what is the expected value of X_n?
    ex_n_odd_formula = "(n^2-1)/2"

    # Question 4: For odd n>30 what is the expected number of times that between
    # the two gifts there are exactly 10 and n-12 friends.
    # This corresponds to a distance of 11 or n-11. The walk is on even integers,
    # so only n-11 is reachable. The formula for the expected number of visits is:
    expected_visits_formula = "4*(n-11)/(n+1)"

    # Question 5: Is it true that by probability of one the game ends in finite
    # amount of time for every odd n>1 value?
    # Yes, for odd n, the random walk of the distance is on a finite, connected
    # state space, so the absorbing state 0 is always reached.
    prob_one_ends = "Yes"

    # Print the answers separated by a comma.
    print(f"{ex_19},{ex_20},{ex_n_odd_formula},{expected_visits_formula},{prob_one_ends}")

solve_and_print()