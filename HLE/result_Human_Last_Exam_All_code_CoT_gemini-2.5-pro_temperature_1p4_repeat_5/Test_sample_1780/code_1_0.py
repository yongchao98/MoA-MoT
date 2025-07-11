import math

def solve_letter_arrangement():
    """
    Calculates the number of ways to arrange the letters "L", "N", "S", and "W"
    based on a specific connection rule.
    """

    # Step 1: Define the properties of the letters.
    # We classify letters by whether their rightmost tip is at the 'Top' or 'Bottom'.
    # - L, N, S: End at the Bottom.
    # - W: Ends at the Top.

    # Step 2: Deduce the connection rule from the problem's examples.
    # - Valid: 'S' (ends Bottom) -> 'W' (starts Top)
    # - Invalid: 'W' (ends Top) -> 'S' (starts Top)
    # The rule is that a letter can only precede another if it ends at the Bottom.

    # Step 3: Apply the rule to find the structure of a valid arrangement.
    # Only 'W' ends at the Top, so it cannot be followed by any other letter.
    # Therefore, 'W' must always be the last letter in the arrangement.

    # Step 4: The problem is now to find the number of ways to arrange the
    # remaining 3 letters (L, N, S) in the first three positions.
    num_letters_to_permute = 3

    # The number of ways to arrange 'n' distinct items is n! (n factorial).
    result = math.factorial(num_letters_to_permute)

    # Print the explanation and the calculation.
    print(f"The letter 'W' must be the last letter in any valid sequence.")
    print(f"The problem reduces to finding the number of permutations of the other {num_letters_to_permute} letters ('L', 'N', 'S').")
    print(f"The number of permutations of {num_letters_to_permute} items is {num_letters_to_permute}!.")
    
    # Building the equation string: 3 * 2 * 1 = 6
    equation_parts = [str(i) for i in range(num_letters_to_permute, 0, -1)]
    equation_str = " * ".join(equation_parts)

    print(f"Calculation: {equation_str} = {result}")

    print("\nTotal number of valid arrangements:")
    print(result)

solve_letter_arrangement()