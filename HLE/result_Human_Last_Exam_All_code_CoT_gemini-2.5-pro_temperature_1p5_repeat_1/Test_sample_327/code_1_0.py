import math

def solve_hackenbush_pieces():
    """
    Calculates the number of pieces for a fraction in Red-Blue Hackenbush
    and displays the corresponding equation.
    """
    # The fraction to represent
    n_orig = 13
    d_orig = 16

    # Use a copy of the numerator for calculations
    n = n_orig
    d = d_orig

    # Store the signed coefficients (1 for Red, -1 for Blue)
    coeffs = []

    # Loop until the number to represent becomes 0
    while n != 0:
        # Multiply the current value by 2
        y_num = 2 * n
        y_den = d
        
        # Determine the color of the next piece
        if y_num > 0:
            # Piece is Red
            sign = 1
            # Subtract 1 from the value y
            # n_new / d_new = y - 1 = (y_num / y_den) - 1 = (y_num - y_den) / y_den
            n = y_num - y_den
            d = y_den
        else: # y_num <= 0
            # Piece is Blue
            sign = -1
            # Add 1 to the value y
            # n_new / d_new = y + 1 = (y_num / y_den) + 1 = (y_num + y_den) / y_den
            n = y_num + y_den
            d = y_den
            
        coeffs.append(sign)
        
        # Simplify the fraction for the next iteration
        if n != 0:
            common_divisor = math.gcd(n, d)
            n //= common_divisor
            d //= common_divisor

    print(f"To represent the number {n_orig}/{d_orig} in red-blue-Hackenbush, we express it as:")
    
    # Build and print the equation string
    equation_parts = []
    for i, s in enumerate(coeffs):
        denominator = 2**(i + 1)
        equation_parts.append(f"({s}) * (1/{denominator})")
    
    print(f"{n_orig}/{d_orig} = " + " + ".join(equation_parts))

    # Determine piece colors
    colors = ["Red" if c == 1 else "Blue" for c in coeffs]
    print(f"This corresponds to a sequence of pieces: {', '.join(colors)}.")

    num_pieces = len(coeffs)
    print(f"\nThe total number of pieces needed is {num_pieces}.")

solve_hackenbush_pieces()