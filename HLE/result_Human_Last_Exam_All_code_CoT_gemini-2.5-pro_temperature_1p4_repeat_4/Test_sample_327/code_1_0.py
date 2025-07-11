def solve_hackenbush_pieces():
    """
    Calculates the number of pieces to represent a fraction in red-blue-Hackenbush.
    
    The method relies on the binary expansion of the fraction. The canonical game
    for a fraction x is the sum of games for its binary components. The number of
    pieces for the game G(1/2^n) is 2n.
    """
    
    numerator = 13
    denominator = 16

    # Denominator must be a power of 2 for this simple construction.
    if (denominator & (denominator - 1)) != 0 or denominator == 0:
        print("Denominator must be a power of 2.")
        return

    # The power of 2 in the denominator determines the length of the binary fraction.
    # e.g., 16 = 2^4, so we need 4 digits.
    power_of_2 = denominator.bit_length() - 1

    # Get the binary representation of the numerator, padded to the required length.
    # e.g., for 13/16, we get "1101"
    binary_representation = bin(numerator)[2:].zfill(power_of_2)

    total_pieces = 0
    calculation_parts = []

    # Iterate through the binary string. Each '1' corresponds to a game G(1/2^n).
    for i, bit in enumerate(binary_representation):
        if bit == '1':
            # n is the position of the bit, starting from 1
            n = i + 1
            
            # The number of pieces for the game G(1/2^n) is 2n
            pieces_for_term = 2 * n
            
            total_pieces += pieces_for_term
            calculation_parts.append(str(pieces_for_term))

    # Print the final equation showing how the total is calculated
    equation = " + ".join(calculation_parts)
    print(f"The total number of pieces is calculated from the components corresponding to the '1's in the binary expansion 0.{binary_representation} of {numerator}/{denominator}.")
    print(f"Pieces = {equation} = {total_pieces}")


solve_hackenbush_pieces()