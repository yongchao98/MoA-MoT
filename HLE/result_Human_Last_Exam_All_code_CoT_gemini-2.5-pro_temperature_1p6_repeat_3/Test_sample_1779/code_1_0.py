def is_language_member(input_string):
    """
    This function implements the decider for the language L where L = {w | |w| is a perfect square}.
    It determines if the length of the input_string is a perfect square.
    The function is guaranteed to halt for any input string.
    """
    n = len(input_string)
    print(f"Testing the input string: '{input_string}'")
    print(f"Step 1: The length of the string is n = {n}.")
    print("------------------------------------------")
    print(f"Step 2: Check if n = {n} is a perfect square.")
    print("We will check integers k = 0, 1, 2, ... and compute k * k.")

    if n < 0:
        # This case is not possible for string lengths but included for completeness.
        print(f"\nResult: Length {n} is negative. It cannot be a perfect square.")
        print(f"Conclusion: The string '{input_string}' is NOT in L.")
        return

    k = 0
    while True:
        # This part models the "final equation" check for each number k.
        square = k * k
        print(f"Checking k = {k}: The equation is {k} * {k} = {square}")

        if square == n:
            print("------------------------------------------")
            print(f"Result: Match found. {n} is a perfect square.")
            print(f"Conclusion: The string '{input_string}' IS in L.")
            break
        elif square > n:
            print("------------------------------------------")
            print(f"Result: Overshot. {square} > {n}. No integer root exists.")
            print(f"Conclusion: The string '{input_string}' is NOT in L.")
            break
        
        k += 1

# --- Example 1: A string whose length (9) is a perfect square ---
is_language_member("aaaaaaaaa")

print("\n=========================================\n")

# --- Example 2: A string whose length (10) is NOT a perfect square ---
is_language_member("bbbbbbbbbb")
