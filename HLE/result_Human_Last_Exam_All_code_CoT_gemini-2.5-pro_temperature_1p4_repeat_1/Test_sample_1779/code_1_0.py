import math

def decider_for_L(w: str):
    """
    This function acts as a decider for the language L.
    It halts on all inputs and determines if the length of the string w
    is a perfect square.

    Args:
        w: The input string.

    Returns:
        True if w is in L, False otherwise.
    """
    k = len(w)
    print(f"--- Checking string \"{w}\" of length {k} ---")

    if k < 0: # Length is always non-negative
        print(f"Length {k} is invalid. String is not in L.")
        return False
        
    # We only need to check integers n up to sqrt(k).
    # If k is a perfect square, its root will be an integer.
    # Otherwise, n*n will become larger than k.
    n = 0
    while True:
        square = n * n
        # This condition demonstrates the core of the decider's logic.
        if square == k:
            print(f"Success: The length {k} is a perfect square.")
            # The prompt requires outputting the numbers in the final equation
            print(f"Equation: {n} * {n} = {k}")
            print(f"Therefore, the string \"{w}\" is in the language L.\n")
            return True
        if square > k:
            print(f"Failure: The length {k} is not a perfect square.")
            print(f"Reached {n} * {n} = {square}, which is greater than {k}.")
            print(f"Therefore, the string \"{w}\" is not in the language L.\n")
            return False
        n += 1

# --- Main execution ---
# Let's test the decider with a few example strings.
# The alphabet can be anything, e.g., {a, b}.
# The decider only cares about the length of the string.
example_strings = [
    "",           # Length 0
    "a",          # Length 1
    "ab",         # Length 2
    "aba",        # Length 3
    "abab",       # Length 4
    "abbababba"   # Length 9
]

for s in example_strings:
    decider_for_L(s)