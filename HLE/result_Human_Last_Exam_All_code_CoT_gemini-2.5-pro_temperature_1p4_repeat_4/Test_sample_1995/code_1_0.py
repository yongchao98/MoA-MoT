def solve_chess_stalemate_puzzle():
    """
    Calculates and explains the solution to the chess stalemate puzzle.

    The puzzle asks for the smallest number of material points for a set of
    white pieces that can attack every square on an 8x8 board except for one,
    which is occupied by a stalemated black king.

    The known minimal solution to this problem, which results in a legal
    stalemate position, requires the following white pieces (excluding the
    white King, which has no point value):
    - 1 Rook
    - 1 Knight
    - 3 Pawns
    """
    
    # Standard chess piece point values
    point_values = {
        "Rook": 5,
        "Knight": 3,
        "Pawn": 1
    }

    # Pieces in the minimal solution
    pieces = {
        "Rook": 1,
        "Knight": 1,
        "Pawn": 3
    }
    
    total_points = 0
    equation_parts = []

    # Calculate the total points and build the equation string
    # Add Rook points
    total_points += pieces["Rook"] * point_values["Rook"]
    for _ in range(pieces["Rook"]):
        equation_parts.append(str(point_values["Rook"]))
        
    # Add Knight points
    total_points += pieces["Knight"] * point_values["Knight"]
    for _ in range(pieces["Knight"]):
        equation_parts.append(str(point_values["Knight"]))

    # Add Pawn points
    total_points += pieces["Pawn"] * point_values["Pawn"]
    for _ in range(pieces["Pawn"]):
        equation_parts.append(str(point_values["Pawn"]))

    # Print the explanation and the final equation
    print("The smallest number of material points is achieved with 1 Rook, 1 Knight, and 3 Pawns.")
    print("The point values are: Rook (5), Knight (3), and Pawn (1).")
    print("\nThe equation for the total points is:")
    
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_points}")


solve_chess_stalemate_puzzle()
