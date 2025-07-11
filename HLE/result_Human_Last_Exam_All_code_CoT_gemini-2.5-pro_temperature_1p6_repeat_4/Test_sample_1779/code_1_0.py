import math

def decide_language_L(w: str):
    """
    This function implements the decider for the language L.
    L = {w | the length of w is a perfect square}.
    The function always halts and decides membership in L.
    """
    n = len(w)

    # An integer n is a perfect square if its square root is an integer.
    # We handle the case for n=0 separately for clarity, though math.isqrt handles it.
    if n < 0:
        # Length cannot be negative, but for robustness we check.
        is_square = False
        sqrt_n = -1 # Placeholder
    elif n == 0:
        is_square = True
        sqrt_n = 0
    else:
        # math.isqrt calculates the integer square root, which is efficient.
        sqrt_n = math.isqrt(n)
        if sqrt_n * sqrt_n == n:
            is_square = True
        else:
            is_square = False

    # Print the result and the reasoning
    if is_square:
        print(f"ACCEPT: The Turing machine T halts on input '{w}'.")
        # The prompt asks to "output each number in the final equation".
        print(f"Reason: Its length is {n}, which is a perfect square ({sqrt_n} * {sqrt_n} = {n}).\n")
    else:
        print(f"REJECT: The Turing machine T does not halt on input '{w}'.")
        print(f"Reason: Its length is {n}, which is not a perfect square.\n")


def main():
    """
    Main function to test the decider with various strings.
    """
    test_strings = [
        "",         # Length 0 (0*0 = 0)
        "a",        # Length 1 (1*1 = 1)
        "xyz",      # Length 3
        "abcd",     # Length 4 (2*2 = 4)
        "computer", # Length 8
        "decidable",# Length 9 (3*3 = 9)
        "turing machine" # Length 14
    ]

    print("Running the decider for language L = {w | |w| is a perfect square}...\n")
    for s in test_strings:
        decide_language_L(s)

if __name__ == "__main__":
    main()