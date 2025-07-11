import math

def decider_for_L(n):
    """
    This function acts as a decider for the language L.
    It determines if a given length 'n' corresponds to a string in L.
    A string w is in L if and only if its length, n, is a perfect square.
    This algorithm is guaranteed to halt for any non-negative integer n.

    Args:
    n (int): The length of an input string.
    """
    print(f"--- Deciding for string length n = {n} ---")

    # The length of a string cannot be negative.
    if n < 0:
        print(f"Error: Length cannot be negative. String is not in L.")
        return

    # We search for an integer k such that k * k = n.
    # We can start our search from k = 0.
    k = 0
    while True:
        # Calculate the square of the current integer k.
        current_square = k * k

        if current_square == n:
            # We found it! n is a perfect square. The string is in L.
            # As requested, output the equation.
            print(f"{k} * {k} = {n}")
            print(f"Conclusion: Length {n} IS a perfect square. The string is in L.")
            return # Halt and accept.

        elif current_square > n:
            # We've gone past n. No integer square root can exist.
            # The string is not in L.
            print(f"Conclusion: Length {n} is NOT a perfect square. The string is not in L.")
            return # Halt and reject.

        # If neither condition is met, check the next integer.
        k = k + 1

# Demonstrate the decider with a list of test lengths.
test_lengths = [9, 12, 16, 25, 30]
for length in test_lengths:
    decider_for_L(length)
