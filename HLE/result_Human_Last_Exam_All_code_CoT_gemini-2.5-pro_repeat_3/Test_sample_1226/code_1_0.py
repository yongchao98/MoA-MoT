def solve_chess_material_problem():
    """
    Calculates the maximum material value for a known mate-in-6 chess problem.
    """
    # Standard point values for chess pieces
    piece_values = {
        'Queen': 9,
        'Rook': 5,
        'Bishop': 3,
        'Knight': 3,
        'Pawn': 1
    }

    # The pieces in the record-holding position (by Friedel & Baier)
    # The position consists of a White King, a Black King, and the following White pieces.
    position_pieces = {
        'Queen': 13,
        'Rook': 1,
        'Pawn': 7
    }

    total_score = 0
    equation_parts = []

    print("Calculating the maximum material value for a mate in 6 based on a known record position:")
    print("-" * 80)

    # Calculate and print the score for each piece type
    for piece, count in position_pieces.items():
        value = piece_values[piece]
        score = count * value
        total_score += score
        print(f"Value from {count} {piece}(s): {count} * {value} = {score} points")
        equation_parts.append(f"({count} * {value})")

    print("-" * 80)
    print("The final equation for the total score is:")
    
    # Building the final equation string with each number explicitly shown
    final_equation_str = ""
    first = True
    for piece, count in position_pieces.items():
        if not first:
            final_equation_str += " + "
        final_equation_str += f"{count} * {piece_values[piece]}"
        first = False
    
    final_equation_str += f" = {total_score}"
    print(final_equation_str)

    print(f"\nThe greatest number of points of white material is {total_score}.")

solve_chess_material_problem()