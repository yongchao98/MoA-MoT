import math

def solve_hackenbush_pieces(num, den):
    """
    Calculates the number of red and blue pieces needed to represent
    a dyadic fraction in red-blue-Hackenbush.
    """
    # Step 1: Check if the denominator is a power of two and find the exponent n.
    # For a fraction num / (2^n), the total number of pieces is n.
    if den <= 0 or (den & (den - 1) != 0):
        print(f"Error: The denominator {den} is not a power of 2.")
        return

    n = int(math.log2(den))
    total_pieces = n

    # Step 2: Convert the numerator to its binary representation, padded to length n.
    # This binary string determines the color of each piece.
    binary_numerator = bin(num)[2:].zfill(n)

    # Step 3: Count the '1's (red pieces) and '0's (blue pieces).
    red_pieces = binary_numerator.count('1')
    blue_pieces = binary_numerator.count('0')

    # Step 4: Print the results and the final equation.
    print(f"The fraction is {num}/{den}.")
    print(f"The binary representation of {num}/{den} is 0.{binary_numerator}.")
    print(f"This representation requires {red_pieces} red pieces (for the '1's) and {blue_pieces} blue pieces (for the '0's).")
    print("\nThe total number of pieces needed is the sum of red and blue pieces.")
    print(f"Total pieces = {red_pieces} + {blue_pieces} = {total_pieces}")


# The specific numbers from the problem.
numerator = 13
denominator = 16

solve_hackenbush_pieces(numerator, denominator)