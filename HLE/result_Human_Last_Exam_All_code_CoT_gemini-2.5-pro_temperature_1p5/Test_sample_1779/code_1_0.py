import math

def perfect_square_decider(w: str):
    """
    This function is an implementation of a decider for the language L where
    L = {w | |w| is a perfect square}.
    It demonstrates that the language is decidable because this algorithm
    always halts and provides a correct yes/no answer.
    """
    n = len(w)
    print(f"Input string: w = '{w}'")
    print(f"Length of w: |w| = {n}")

    # The core of the decision algorithm: check if n is a perfect square.
    # This check is finite and guaranteed to terminate.
    if n < 0:
        # This case is not possible for string lengths but included for completeness.
        is_perfect_square = False
    elif n == 0:
        is_perfect_square = True
        root = 0
    else:
        # A simple and efficient way to check for a perfect square.
        # This is equivalent to our iterative search algorithm.
        root = int(math.sqrt(n))
        if root * root == n:
            is_perfect_square = True
        else:
            is_perfect_square = False

    # Output the final decision.
    if is_perfect_square:
        print(f"Result: Yes, the length {n} is a perfect square.")
        # As requested, output the numbers in the final equation.
        print(f"Equation: {root} * {root} = {n}")
    else:
        print(f"Result: No, the length {n} is not a perfect square.")

if __name__ == '__main__':
    # Test with a string whose length is a perfect square.
    perfect_square_decider("HelloWorld") # Length is 10, not a perfect square.
    print("\n" + "="*30 + "\n")

    # Test with a string whose length is a perfect square.
    perfect_square_decider("123456789") # Length is 9, a perfect square.
    print("\n" + "="*30 + "\n")

    # Test with the empty string.
    perfect_square_decider("") # Length is 0, a perfect square.
    print("\n" + "="*30 + "\n")

    # Test with a string of length 1.
    perfect_square_decider("x") # Length is 1, a perfect square.