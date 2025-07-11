import collections

def solve_chess_puzzle():
    """
    Calculates and prints the solution for the Diagonal Corridor Mate puzzle.
    """
    # Define standard piece values.
    piece_values = {
        "Pawn": 1,
        "Knight": 3,
        "Bishop": 3,
        "Rook": 5,
        "Queen": 9
    }

    # To create a true "Diagonal Corridor Mate", the Black King on h8 must be
    # trapped by a corridor of its own pieces.
    # The minimum number of pieces to form a corridor is two.

    # 1. Determine the White checking piece.
    # A piece on the a1-h8 diagonal is needed. A Bishop is the minimal choice.
    white_piece_name = "Bishop"
    white_piece_value = piece_values[white_piece_name]
    white_pieces_list = [f"White {white_piece_name}"]

    # 2. Determine the Black corridor pieces.
    # The White Bishop attacks g7. The king's other escape squares, h7 and g8,
    # must be blocked by its own pieces. Two Pawns are the minimal choice.
    black_piece_names = ["Pawn", "Pawn"]
    black_piece_values = [piece_values[name] for name in black_piece_names]
    
    # Sort black pieces by value to ensure correct order if pieces were different.
    # In this case, they are the same.
    black_pieces_list = [f"Black {name}" for name in sorted(black_piece_names, key=lambda p: piece_values[p])]

    # 3. Combine lists for the final answer.
    final_piece_list = white_pieces_list + black_pieces_list
    final_answer_string = ", ".join(final_piece_list)

    # 4. Calculate the total value for the equation.
    total_value = white_piece_value + sum(black_piece_values)

    # 5. Output the explanation, equation, and final answer.
    print("The minimal configuration for a 'Diagonal Corridor Mate' requires 3 additional pieces.")
    print("This ensures the king is trapped by a 'corridor' of its own pieces, as the pattern name implies.")
    print("\n- A White Bishop delivers the check along the diagonal and guards one escape square (g7).")
    print("- Two Black Pawns (on h7 and g8) form the corridor, blocking the other escape squares.")
    
    # As requested, output each number in the final equation for total piece value.
    value_equation_str = f"{white_piece_value} + {' + '.join(map(str, black_piece_values))} = {total_value}"
    print("\nThe equation for the minimum total piece value is:")
    print(value_equation_str)

    print("\nThe final comma-separated list of additional pieces is:")
    print(final_answer_string)

solve_chess_puzzle()