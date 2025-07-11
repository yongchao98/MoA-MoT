import collections

def solve_chess_puzzle():
    """
    This function determines and prints the minimal piece set for a specific
    Diagonal Corridor Mate.
    """

    # Piece values for sorting purposes.
    piece_values = {
        "Pawn": 1,
        "Knight": 3,
        "Bishop": 3,
        "Rook": 5,
        "Queen": 9
    }

    # The minimal set of additional pieces determined by the logic.
    # White pieces deliver the mate.
    # Black pieces form the corridor, trapping their own king.
    white_pieces_solution = ["White Bishop"]
    black_pieces_solution = ["Black Pawn", "Black Pawn"]

    # Function to get the base piece name for sorting.
    def get_value(piece_name_with_color):
        base_name = piece_name_with_color.split()[-1]
        return piece_values.get(base_name, 0)

    # Sort each list by piece value in ascending order.
    # This is for correctness, although the lists are already simple.
    white_pieces_solution.sort(key=get_value)
    black_pieces_solution.sort(key=get_value)

    # Combine the lists as per the output requirement (White, then Black).
    final_pieces_list = white_pieces_solution + black_pieces_solution

    # Print the final comma-separated string.
    print(", ".join(final_pieces_list))

solve_chess_puzzle()