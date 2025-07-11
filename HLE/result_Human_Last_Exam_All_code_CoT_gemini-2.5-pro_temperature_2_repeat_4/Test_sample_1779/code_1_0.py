import math

def check_perfect_square(n):
    """
    This function implements an algorithm to decide if a non-negative integer n
    is a perfect square. This demonstrates the decidability of the language L,
    where L = {w | length(w) is a perfect square}.

    The algorithm always halts, proving L is decidable.

    Args:
        n: A non-negative integer representing the length of an input string.
    """
    if not isinstance(n, int) or n < 0:
        print(f"Error: Input must be a non-negative integer. Received: {n}")
        return

    # A highly efficient way to check for a perfect square.
    # We find the integer square root and check if its square equals n.
    # This is equivalent to the Turing machine algorithm of checking i*i
    # but much faster in Python.
    root = math.isqrt(n)

    # The "final equation" part as requested by the prompt.
    if root * root == n:
        print(f"HALT and ACCEPT: The length {n} is a perfect square. Equation: {root} * {root} = {n}")
    else:
        # To show why it's not a perfect square, we find the two perfect squares it lies between.
        lower_root = root
        upper_root = root + 1
        lower_square = lower_root * lower_root
        upper_square = upper_root * upper_root
        print(f"HALT and REJECT: The length {n} is NOT a perfect square. It lies between {lower_root}*{lower_root}={lower_square} and {upper_root}*{upper_root}={upper_square}.")

if __name__ == '__main__':
    # Demonstrate the decision procedure for lengths from 0 to 50.
    print("Testing if a Turing machine T halts on strings of length n (0 to 50)...")
    print("T halts if and only if length(w) is a perfect square.\n")
    for length in range(51):
        check_perfect_square(length)
