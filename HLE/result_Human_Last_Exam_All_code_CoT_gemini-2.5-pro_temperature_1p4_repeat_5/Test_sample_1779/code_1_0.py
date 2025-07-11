import math

def turing_machine_decider(w: str):
    """
    This function simulates a decider for the language L where a string's
    length must be a perfect square.
    It demonstrates that the language is decidable by providing a clear,
    always-halting algorithm.
    """
    print(f"Analyzing input string w: '{w}'")
    
    # Step 1: Calculate the length of the string w.
    n = len(w)
    print(f"The length of w is |w| = {n}.")

    # Step 2: Determine if n is a perfect square.
    print(f"Checking if n = {n} is a perfect square.")
    
    # An effective way to check for a perfect square is to find the
    # integer square root and see if squaring it gives back the original number.
    if n < 0:
        is_square = False
        root = -1 # Placeholder
    else:
        root = int(math.sqrt(n))
        is_square = (root * root == n)

    # Step 3: Conclude and show the final check/equation.
    print(f"The integer square root of {n} is {root}.")
    print(f"Final equation check: {root} * {root} = {root * root}")

    if is_square:
        print(f"\nResult: The length {n} is a perfect square.")
        print(f"Therefore, the string '{w}' is in the language L, and our decider ACCEPTS it.")
    else:
        print(f"\nResult: The length {n} is not a perfect square.")
        print(f"Therefore, the string '{w}' is not in the language L, and our decider REJECTS it.")


# --- Main Execution ---
# We will use an example string with a length of 9, which is a perfect square.
input_string = "perfect!!"
turing_machine_decider(input_string)
