import math

def solve_hackenbush_pieces():
    """
    Calculates the number of pieces needed to represent a dyadic rational
    in red-blue-Hackenbush.
    """
    # The given number is 13/16.
    num = 13
    den = 16

    # 1. Calculate the number of pieces for the integer part.
    # This is the absolute value of the integer part of the number.
    integer_part = num // den
    integer_pieces = abs(integer_part)

    # 2. Calculate the number of pieces for the fractional part.
    # For a fraction m/(2^n), this is 'n'. 'n' is the base-2 log of the denominator.
    # We first verify the denominator is a power of two.
    if not (den > 0 and (den & (den - 1) == 0)):
        print(f"Error: The denominator {den} is not a power of 2.")
        return

    fractional_pieces = int(math.log2(den))

    # 3. The total number of pieces is the sum of the two parts.
    total_pieces = integer_pieces + fractional_pieces
    
    # Print the final equation showing how the total is calculated.
    print(f"To represent the number {num}/{den}:")
    print(f"Pieces for integer part ({integer_part}) + Pieces for fractional part (from denominator {den}) = Total Pieces")
    print(f"{integer_pieces} + {fractional_pieces} = {total_pieces}")

solve_hackenbush_pieces()