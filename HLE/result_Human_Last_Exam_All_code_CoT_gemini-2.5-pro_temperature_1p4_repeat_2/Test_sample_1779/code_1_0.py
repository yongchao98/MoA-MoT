import math

def decide_language_L(w: str):
    """
    This function acts as a decider for the language L.

    The problem states:
    1. A Turing machine T halts on input w if and only if the length of w is a perfect square.
    2. A language L contains a string w if and only if T halts on w.

    Therefore, a string w is in L if and only if its length is a perfect square.
    This function determines if a given string w belongs to L by checking
    if its length is a perfect square. This check is an always-halting
    procedure, which proves that the language L is decidable.
    """
    n = len(w)

    print(f"Input string: '{w}'")
    print(f"Length of the string is: {n}")

    # The core of the decider: check if n is a perfect square.
    if n < 0:
        # This case is not possible for string length.
        is_square = False
    else:
        root = int(math.sqrt(n))
        is_square = (root * root == n)

    print(f"Checking if the length ({n}) is a perfect square...")
    if is_square:
        # Fulfilling the request to show the equation
        print(f"Result: Yes, because {root} * {root} = {n}.")
        print(f"Since the length is a perfect square, T halts on '{w}'.")
        print(f"Therefore, '{w}' IS in the language L.")
    else:
        # Fulfilling the request to show numbers for the non-square case
        lower_root = int(math.sqrt(n))
        upper_root = lower_root + 1
        print(f"Result: No. The number {n} is not a perfect square.")
        print(f"It falls between the two consecutive squares: {lower_root}*{lower_root}={lower_root**2} and {upper_root}*{upper_root}={upper_root**2}.")
        print(f"Since the length is not a perfect square, T does not halt on '{w}'.")
        print(f"Therefore, '{w}' IS NOT in the language L.")
    print("-" * 40)


# --- Main execution with examples ---

# Example 1: A string with a length that is a perfect square (9).
decide_language_L("aaaaaaaaa")

# Example 2: A string with a length that is not a perfect square (10).
decide_language_L("bbbbbbbbbb")

# Example 3: The empty string (length 0, which is a perfect square).
decide_language_L("")

# Example 4: A string of length 1 (which is a perfect square).
decide_language_L("x")