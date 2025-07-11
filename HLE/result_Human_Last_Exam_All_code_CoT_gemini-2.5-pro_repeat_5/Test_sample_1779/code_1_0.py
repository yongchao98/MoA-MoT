import math

def decider_for_L(input_string):
    """
    This function acts as a decider for the language L.
    L is the set of strings whose length is a perfect square.
    The function determines if the input string's length is a perfect square.
    It always halts and prints a decision.
    """
    k = len(input_string)
    print(f"Analyzing input string 'w' of length k = {k}.")

    # A simple and efficient way to check for a perfect square.
    # This avoids a loop but follows the same decidable logic.
    if k < 0:
        is_perfect = False
    elif k == 0:
        is_perfect = True
        root = 0
    else:
        root = int(math.sqrt(k))
        if root * root == k:
            is_perfect = True
        else:
            is_perfect = False

    # Output the result based on the decision.
    if is_perfect:
        print(f"The length {k} is a perfect square.")
        # Output the equation as requested.
        print(f"The final equation is: {root} * {root} = {k}")
        print("Result: The Turing machine T halts on w, so w IS in language L.")
    else:
        print(f"The length {k} is not a perfect square.")
        # If not a perfect square, show why by checking the nearest integer root.
        print(f"Checking the integer part of the square root: {root} * {root} = {root*root}.")
        print(f"And the next integer's square: {root+1} * {root+1} = {(root+1)*(root+1)}.")
        print("Result: The Turing machine T does not halt on w, so w IS NOT in language L.")

if __name__ == '__main__':
    # --- Test Case 1: A string whose length IS a perfect square ---
    # The string can contain any character; only its length matters.
    w1 = "x" * 49
    print("--- Test Case 1 ---")
    decider_for_L(w1)

    print("\n" + "="*50 + "\n")

    # --- Test Case 2: A string whose length IS NOT a perfect square ---
    w2 = "y" * 50
    print("--- Test Case 2 ---")
    decider_for_L(w2)