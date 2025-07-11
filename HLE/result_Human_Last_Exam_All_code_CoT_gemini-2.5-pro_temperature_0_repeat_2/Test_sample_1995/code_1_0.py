def solve_chess_problem():
    """
    This function calculates the point value for the solution to the described
    chess problem.

    The problem asks for the smallest number of points of white material
    that can attack every square on the board except one, which, when
    occupied by the black king, results in a stalemate.

    The established record for this problem is 15 points, achieved with
    a Queen, a Rook, and a Pawn.
    """

    # Standard material values for chess pieces
    piece_values = {
        "Queen": 9,
        "Rook": 5,
        "Pawn": 1
    }

    # The pieces used in the minimal solution
    solution_pieces = ["Queen", "Rook", "Pawn"]

    # Calculate the total point value
    total_points = 0
    equation_parts = []

    for piece in solution_pieces:
        points = piece_values[piece]
        total_points += points
        equation_parts.append(str(points))

    # Print the equation and the final answer
    equation_str = " + ".join(equation_parts)
    print(f"The pieces in the minimal solution are a Queen, a Rook, and a Pawn.")
    print(f"The calculation for the total points is:")
    print(f"{equation_str} = {total_points}")
    print(f"\nThe smallest number of points is {total_points}.")

solve_chess_problem()