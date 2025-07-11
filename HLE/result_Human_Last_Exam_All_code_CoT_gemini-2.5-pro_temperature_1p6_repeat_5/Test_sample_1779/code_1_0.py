import math

def decider_for_L(input_string):
    """
    This function simulates a decider for the language L.
    L is the set of strings whose length is a perfect square.
    The decider must always halt and provide a 'yes' or 'no' answer.
    """
    n = len(input_string)
    print(f"Analyzing input string: '{input_string}'")
    print(f"String length is: {n}")

    # Handle the trivial case of a negative length, which is impossible for strings
    # but good practice for number theory functions.
    if n < 0:
        print("Result: String length cannot be negative. String is not in L.")
        return

    # A simple and always-halting way to check for a perfect square.
    # We find the integer closest to the square root and check if its square is n.
    # This avoids loops for efficiency but demonstrates the same decidable principle.
    sqrt_n = int(math.sqrt(n))

    # This is the 'final equation' check.
    # If sqrt_n * sqrt_n == n, then n is a perfect square.
    if sqrt_n * sqrt_n == n:
        print(f"Check: The number {n} is a perfect square.")
        print(f"Final Equation: {sqrt_n} * {sqrt_n} = {n}")
        print("Result: The Turing machine T halts on this input. The string is in language L.")
    else:
        print(f"Check: The number {n} is not a perfect square.")
        print("Result: The Turing machine T does not halt on this input. The string is not in language L.")
    print("-" * 20)

# --- Main Execution ---
# Let's test our decider with a few example strings.

# Test Case 1: Length 9 (a perfect square)
decider_for_L("#########")

# Test Case 2: Length 7 (not a perfect square)
decider_for_L("abcdefg")

# Test Case 3: Length 0 (a perfect square, 0*0=0)
decider_for_L("")

# Test Case 4: Length 16 (a perfect square)
decider_for_L("1234567890123456")
