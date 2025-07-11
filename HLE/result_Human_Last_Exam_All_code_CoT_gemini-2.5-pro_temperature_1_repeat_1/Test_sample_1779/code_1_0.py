def check_perfect_square_language(input_string):
    """
    This function simulates a decider for the language L where
    w is in L if and only if len(w) is a perfect square.
    The function always halts and decides membership in L.
    """
    n = len(input_string)
    print(f"Input string: '{input_string}'")
    print(f"Step 1: Calculate the length of the string.")
    print(f"Length n = {n}\n")

    print(f"Step 2: Check if n = {n} is a perfect square.")
    print("We will check k*k for k = 0, 1, 2, ...")

    k = 0
    while True:
        k_squared = k * k
        print(f"Checking k = {k}: The equation is {k} * {k} = {k_squared}")

        if k_squared == n:
            print(f"\nResult: Since {k_squared} == {n}, the length is a perfect square.")
            print(f"The string '{input_string}' is in the language L.")
            break
        elif k_squared > n:
            print(f"\nResult: Since {k_squared} > {n}, we have passed {n}.")
            print("The length is not a perfect square.")
            print(f"The string '{input_string}' is NOT in the language L.")
            break
        k += 1

    print("\nBecause this procedure is guaranteed to halt for any input, the language L is decidable.")

# Example usage with a string of length 49
example_string = "a" * 49
check_perfect_square_language(example_string)

print("-" * 20)

# Example usage with a string of length 50
example_string_2 = "b" * 50
check_perfect_square_language(example_string_2)
