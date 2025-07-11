from fractions import Fraction

def solve_hackenbush_number(n, d):
    """
    Calculates the number of pieces and the equation for a dyadic rational
    in Red-Blue Hackenbush.
    
    Args:
        n (int): The numerator of the fraction.
        d (int): The denominator of the fraction (must be a power of 2).
    """
    if d & (d - 1) != 0 and d != 0:
        print("Error: The denominator must be a power of 2.")
        return

    target_fraction = Fraction(n, d)
    current_value = target_fraction

    # s_values will store +1 for Red pieces and -1 for Blue pieces.
    s_values = []
    
    while current_value != 0:
        if current_value > 0:
            s_values.append(1)
            current_value = 2 * current_value - 1
        else: # current_value < 0
            s_values.append(-1)
            current_value = 2 * current_value + 1

    num_pieces = len(s_values)

    print(f"The number {target_fraction} can be represented by {num_pieces} pieces in Red-Blue Hackenbush.")
    print("The corresponding equation is:")
    
    # Construct the list of elements to be printed for the equation
    equation_elements = []
    for i, s in enumerate(s_values):
        term_denom = 2**(i + 1)
        if i == 0:
            # Handle the first term, which doesn't have a preceding operator unless negative
            if s < 0:
                equation_elements.append("-")
            equation_elements.append(f"1/{term_denom}")
        else:
            # Handle subsequent terms with their operators
            if s > 0:
                equation_elements.append("+")
            else:
                equation_elements.append("-")
            equation_elements.append(f"1/{term_denom}")
            
    equation_elements.append("=")
    equation_elements.append(str(target_fraction))
    
    # Print the equation with spaces between each element
    print(*equation_elements)

# Solve for the given number 13/16
solve_hackenbush_number(13, 16)