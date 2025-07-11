def hackenbush_pieces(numerator, denominator):
    """
    Calculates the number of red and blue pieces needed to represent a
    dyadic rational number in Red-Blue-Hackenbush.
    """
    if not (denominator > 0 and denominator & (denominator - 1) == 0):
        print("Error: The denominator must be a power of 2.")
        return

    num = numerator
    den = denominator
    
    red_pieces = 0
    blue_pieces = 0
    binary_representation = "0."
    
    # The loop continues until the fraction is resolved.
    # For a dyadic rational m / 2^n, this will take n steps.
    while num > 0:
        num *= 2
        if num >= den:
            binary_representation += "1"
            red_pieces += 1
            num -= den
        else:
            binary_representation += "0"
            blue_pieces += 1
            
    total_pieces = red_pieces + blue_pieces
    
    print(f"To represent the fraction {numerator}/{denominator}, we use its binary form: {binary_representation}")
    print(f"Number of red pieces (1s): {red_pieces}")
    print(f"Number of blue pieces (0s): {blue_pieces}")
    # Final equation showing the sum
    print(f"Total number of pieces: {red_pieces} + {blue_pieces} = {total_pieces}")

# Input numbers for the problem
numerator = 13
denominator = 16

hackenbush_pieces(numerator, denominator)