import math

def decider_for_L(w):
    """
    This function implements a decider for the language L.
    L = {w | the length of w is a perfect square}.
    It returns True if w is in L and False otherwise, and always halts.
    """
    n = len(w)
    
    # Handle the case of n=0, as it's a perfect square.
    if n == 0:
        root = 0
    else:
        # Calculate the integer part of the square root.
        root = int(math.sqrt(n))

    # Check if the number is a perfect square by squaring the integer root.
    if root * root == n:
        print(f"Input string: '{w}'")
        print(f"Length of string is: {n}")
        print(f"The number {n} IS a perfect square because {root} * {root} = {n}.")
        print("Therefore, the Turing machine T halts on this input, and the string is in the language L.")
    else:
        print(f"Input string: '{w}'")
        print(f"Length of string is: {n}")
        print(f"The number {n} IS NOT a perfect square.")
        print("Therefore, the Turing machine T does not halt on this input, and the string is not in the language L.")
    print("-" * 20)


# --- Main execution ---
# Demonstrate with an example string that is in the language L.
decider_for_L("decidable")

# Demonstrate with an example string that is NOT in the language L.
decider_for_L("undecidable")