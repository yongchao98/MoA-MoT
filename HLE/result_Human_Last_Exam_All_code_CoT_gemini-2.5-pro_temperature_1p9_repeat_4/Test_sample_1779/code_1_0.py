def decide_perfect_square_length(w: str):
    """
    This function simulates a decider for the language L where a string's
    length must be a perfect square. It demonstrates that the language is decidable
    by implementing an algorithm that always halts.
    """
    n = len(w)
    print(f"Analyzing input string: \"{w}\"")
    print(f"Length of the string is n = {n}")
    print("---------------------------------")
    print("Starting decision process...")

    # A decider must halt. This algorithm checks squares sequentially.
    i = 0
    while True:
        square = i * i
        print(f"Checking if n matches i*i, where i={i}: {i}*{i} = {square}")

        if square == n:
            # The condition is met. The string is in the language.
            print("\nResult: HALT and ACCEPT.")
            print(f"The length {n} is a perfect square.")
            print(f"Final Equation: {n} = {i} * {i}")
            return True

        if square > n:
            # We have surpassed n. It can't be a perfect square.
            # The string is not in the language.
            print(f"\nResult: HALT and REJECT.")
            print(f"The check {square} > {n} is true, so {n} is not a perfect square.")
            return False

        # If neither condition is met, continue to the next integer.
        i += 1

# --- Main Execution ---

# Test Case 1: A string whose length (9) is a perfect square.
print("--- Test Case 1 ---")
decide_perfect_square_length("decidable")
print("\n=================================\n")


# Test Case 2: A string whose length (12) is not a perfect square.
print("--- Test Case 2 ---")
decide_perfect_square_length("Turing test")
print("\n=================================\n")

# Test Case 3: A string whose length (0) is a perfect square.
print("--- Test Case 3 ---")
decide_perfect_square_length("")