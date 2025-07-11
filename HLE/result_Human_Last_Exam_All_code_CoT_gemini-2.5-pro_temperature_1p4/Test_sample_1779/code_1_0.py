import math

def decider_for_L(input_string):
    """
    This function acts as a decider for the language L.
    It determines if a string's length is a perfect square.
    It always halts and provides a clear 'yes' or 'no' answer.
    """
    n = len(input_string)
    print(f"Input string: '{input_string}'")
    print(f"String length (n): {n}")

    # A simple case for non-negative lengths
    if n < 0:
        print(f"Result: The length {n} is negative. String is not in L.")
        return

    # Start checking for an integer square root from k=0
    k = 0
    while True:
        square = k * k
        # Check if k*k is the target length n
        if square == n:
            print(f"Result: T halts on this input. The string is in L.")
            # The final equation as requested
            print(f"Because its length is a perfect square: {k} * {k} = {n}")
            return
        # If k*k exceeds n, n cannot be a perfect square
        elif square > n:
            print(f"Result: T does not halt on this input. The string is not in L.")
            print(f"Because its length {n} is not a perfect square.")
            return
        # Continue to the next integer
        k += 1

# --- Main execution ---
# Demonstrate the decider with a few examples

# Example 1: A string whose length is a perfect square (9)
w1 = "110101001"
decider_for_L(w1)
print("-" * 30)

# Example 2: A string whose length is not a perfect square (12)
w2 = "a" * 12
decider_for_L(w2)
print("-" * 30)

# Example 3: A string of length 0 (which is a perfect square)
w3 = ""
decider_for_L(w3)
print("-" * 30)

# Example 4: A string of length 1 (which is a perfect square)
w4 = "x"
decider_for_L(w4)
