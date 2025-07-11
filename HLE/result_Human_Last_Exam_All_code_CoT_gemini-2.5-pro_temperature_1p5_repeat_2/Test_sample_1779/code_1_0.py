import math

def decider_for_L(w):
    """
    This function simulates a decider for the language L, where L is the set
    of all strings whose length is a perfect square.

    The decider must always halt and output either ACCEPT or REJECT.
    """
    # Step 1: Compute the length of the input string w.
    n = len(w)
    print(f"Input string: '{w}'")
    print(f"Length of string: n = {n}")

    # Step 2 & 3: Check if n is a perfect square. This is a finite process
    # that is guaranteed to halt.
    if n < 0: # Length of a string cannot be negative, but for completeness.
        is_square = False
        root = -1
    elif n == 0:
        is_square = True
        root = 0
    else:
        # A simple and efficient way to check for a perfect square.
        root = int(math.sqrt(n))
        is_square = (root * root == n)

    # Output the result based on the check.
    if is_square:
        print(f"The length {n} is a perfect square.")
        # Output the equation as requested.
        print(f"Equation: {root} * {root} = {n}")
        print("Result: ACCEPT")
    else:
        print(f"The length {n} is not a perfect square.")
        print("Result: REJECT")
    
    # Separator for clarity
    print("-" * 20)

def main():
    """
    Demonstrates the decider with a few example strings.
    """
    print("Simulating a decider for the language L = {w | |w| is a perfect square}.\n")

    # A string whose length (9) is a perfect square.
    decider_for_L("algorithm")

    # A string whose length (12) is not a perfect square.
    decider_for_L("computations")
    
    # An empty string, whose length (0) is a perfect square (0*0=0).
    decider_for_L("")

    # A string of length 4, which is a perfect square.
    decider_for_L("bits")

if __name__ == "__main__":
    main()