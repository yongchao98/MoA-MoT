def decide_language_for_perfect_square_length(w):
    """
    This function decides if a string 'w' belongs to the language L,
    where L = {w | the length of w is a perfect square}.
    It simulates the decider Turing machine M described above.
    """
    n = len(w)
    print(f"Input string: w = '{w}'")
    print(f"Length of w: |w| = {n}")
    print("---------------------------------")
    print(f"Checking if {n} is a perfect square...")

    # The loop simulates the systematic check performed by the Turing machine.
    # It starts checking integers k = 0, 1, 2, ...
    k = 0
    while True:
        square = k * k
        print(f"Checking k = {k}: Is {k} * {k} == {n}? Result: {square}")
        
        # This is the "accept" condition.
        if square == n:
            print(f"\nConclusion: YES, {n} is a perfect square.")
            print(f"The final equation is {k} * {k} = {n}.")
            print(f"Therefore, the string '{w}' is in the language L.")
            return # The decider halts and accepts.
        
        # This is the "reject" condition.
        if square > n:
            print(f"\nConclusion: NO, {n} is not a perfect square because {square} > {n}.")
            print(f"Therefore, the string '{w}' is NOT in the language L.")
            return # The decider halts and rejects.
        
        # If neither condition is met, continue to the next integer.
        k += 1

# --- Example 1: A string whose length is a perfect square (length = 9) ---
# A string in L. The decider should halt and accept.
print("Running Example 1:")
decide_language_for_perfect_square_length("PythonFun")
print("\n================================================\n")


# --- Example 2: A string whose length is NOT a perfect square (length = 12) ---
# A string not in L. The decider should halt and reject.
print("Running Example 2:")
decide_language_for_perfect_square_length("TuringMachine")
