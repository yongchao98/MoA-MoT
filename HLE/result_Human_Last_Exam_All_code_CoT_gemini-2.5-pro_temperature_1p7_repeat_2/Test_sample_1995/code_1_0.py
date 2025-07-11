def solve_chess_stalemate_problem():
    """
    Calculates and displays the material value for the optimal solution
    to the described chess problem.
    """
    # Standard material values in chess (King has no point value)
    piece_values = {
        'Queen': 9,
        'Rook': 5,
        'Bishop': 3,
        'Knight': 3,
        'Pawn': 1
    }

    # Pieces in the optimal solution (excluding the King)
    white_pieces = {
        'Queen': 1,
        'Rook': 1,
        'Bishop': 2,
        'Pawn': 2
    }

    print("The solution uses the following white pieces (excluding the King):")
    for piece, count in white_pieces.items():
        print(f"- {count} x {piece} ({piece_values[piece]} points each)")
    
    print("\nThe calculation for the total material value is:")
    
    total_points = 0
    calculation_parts = []
    
    # Build the equation string step-by-step
    for piece, count in white_pieces.items():
        value = piece_values[piece]
        for _ in range(count):
            total_points += value
            calculation_parts.append(str(value))

    # Print the equation with each number
    equation = " + ".join(calculation_parts)
    print(f"{equation} = {total_points}")

    print(f"\nThe smallest number of points of white material is {total_points}.")

solve_chess_stalemate_problem()
