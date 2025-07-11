import math

def is_language_member(w: str):
    """
    This function acts as a decider for the language L, where L is the set of
    strings whose length is a perfect square. It demonstrates that the language
    is decidable because this algorithm always halts and gives a clear answer.
    
    Args:
        w: The input string.
    """
    length = len(w)
    
    print(f"Input string w: '{w}'")
    print(f"Step 1: Calculate the length of w.")
    print(f"Length n = {length}")
    print("-" * 20)
    print(f"Step 2: Check if n = {length} is a perfect square.")
    
    # An efficient way to check for a perfect square.
    # This check is guaranteed to complete in a finite number of steps.
    if length < 0:
        # The length of a string is always non-negative.
        is_perfect_square = False
    elif length == 0:
        is_perfect_square = True
        root = 0
    else:
        # Calculate the integer part of the square root.
        root = int(math.sqrt(length))
        # Check if squaring the root gives back the original length.
        if root * root == length:
            is_perfect_square = True
        else:
            is_perfect_square = False

    print("Step 3: Conclude and halt.")
    if is_perfect_square:
        print(f"The final equation is: {length} = {root} * {root}")
        print(f"Result: The length {length} is a perfect square.")
        print("Therefore, the string w is in the language L, and the Turing machine T halts on w.")
    else:
        # Show why it's not a perfect square
        print(f"The check fails: {root} * {root} = {root * root}, which does not equal {length}.")
        print(f"Result: The length {length} is not a perfect square.")
        print("Therefore, the string w is not in the language L, and the Turing machine T does not halt on w.")

# --- Main execution ---
# Let's test with a string whose length is 16 (a perfect square).
test_string = "a" * 16
is_language_member(test_string)
