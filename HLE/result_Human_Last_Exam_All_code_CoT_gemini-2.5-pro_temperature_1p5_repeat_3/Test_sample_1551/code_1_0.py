def solve_puzzle():
    """
    This function calculates and prints the answers to the puzzle.
    The answers are derived from modeling the gift-passing game as a random walk.
    """

    # Answer 1: Expected value of X_19
    # For odd n, the expected time is (n^2-1)/2.
    n_19 = 19
    # The numbers in the equation are 19, 2, 1. The result is 180.
    # Calculation: (19*19 - 1) / 2 = (361 - 1) / 2 = 360 / 2 = 180
    expected_time_19 = (n_19**2 - 1) // 2

    # Answer 2: Expected value of X_20
    # For even n, the game never ends.
    expected_time_20 = "infinity"

    # Answer 3: Expected value of X_n for odd n > 1
    # The general formula is (n^2-1)/2. The numbers are n, 2, 1.
    general_expected_time_odd_n = "(n^2-1)/2"

    # Answer 4: For odd n > 30, expected number of times the distance is 11 (10 friends between).
    # The formula is 2(n-13). The numbers are 2, n, 13.
    expected_visits = "2(n-13)"

    # Answer 5: Is it true that the game ends in finite time with probability 1 for every odd n > 1?
    # Yes, because it's a random walk on a finite state space with a reachable absorbing state.
    finite_time_prob_one = "Yes"

    # Combine all answers into a single string, separated by commas.
    final_answer = f"{expected_time_19},{expected_time_20},{general_expected_time_odd_n},{expected_visits},{finite_time_prob_one}"

    print(final_answer)

solve_puzzle()