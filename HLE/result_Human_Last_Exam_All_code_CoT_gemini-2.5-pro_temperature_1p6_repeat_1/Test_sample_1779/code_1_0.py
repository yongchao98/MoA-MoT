import math

def decider_for_L(w: str):
    """
    This function acts as a decider for the language L = {w | |w| is a perfect square}.
    It determines if the Turing machine T (which halts on strings of perfect square length)
    would halt on the input string w.

    Args:
        w: The input string.
    """
    n = len(w)
    print(f"Input string: '{w}'")
    print(f"Length of string, n = {n}")

    # The algorithm to decide if n is a perfect square.
    # This algorithm is guaranteed to halt for any non-negative integer n.
    if n < 0:
        # This case is not possible for string lengths but included for completeness.
        is_square = False
    elif n == 0:
        is_square = True
        root = 0
    else:
        # Calculate the integer part of the square root
        root = int(math.sqrt(n))
        # Check if the square of the integer root equals the original number
        if root * root == n:
            is_square = True
        else:
            is_square = False

    # Output the result of the decision procedure.
    if is_square:
        # Show the "equation" demonstrating it's a perfect square
        print(f"Decision: Length {n} IS a perfect square because {root} * {root} = {n}.")
        print("Result: Therefore, the language L contains this string, and T halts.\n")
    else:
        # Show why it's not a perfect square
        print(f"Decision: Length {n} IS NOT a perfect square.")
        print(f"Result: Therefore, the language L does not contain this string, and T does not halt.\n")


if __name__ == "__main__":
    # A list of test strings to check
    test_strings = [
        "ninechars", # length 9, a perfect square
        "a",         # length 1, a perfect square
        "hello",     # length 5, not a perfect square
        "",          # length 0, a perfect square
        "decidable", # length 9, a perfect square
        "computer"   # length 8, not a perfect square
    ]

    for s in test_strings:
        decider_for_L(s)