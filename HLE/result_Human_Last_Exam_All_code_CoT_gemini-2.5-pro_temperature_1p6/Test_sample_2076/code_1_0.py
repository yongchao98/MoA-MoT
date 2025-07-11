def solve_chess_puzzle():
    """
    This function determines the minimum pieces for a Diagonal Corridor Mate.
    
    The problem specifies:
    - White King on a1, Black King on h8.
    - White delivers checkmate.
    - The pattern is a Diagonal Corridor Mate.
    - The goal is the minimum number of additional pieces.
    - If multiple solutions have the same number of pieces, choose the one with minimum total piece value.

    Reasoning:
    1. A White Bishop on g7 checks the Black King on h8 and controls its escape squares g8 and g7.
    2. The Bishop on g7 must be protected from being captured by the King.
    3. The cheapest piece to protect the Bishop on g7 is a White Pawn on f6.
    4. This 2-piece solution (White Pawn, White Bishop) is minimal in piece count and value.

    The required output format is a comma-separated string, listing White pieces first,
    then Black, with each sub-list sorted by piece value.
    - White pieces: Pawn (value 1), Bishop (value 3). Sorted: "White Pawn", "White Bishop".
    - Black pieces: None.
    """

    # The additional pieces required for the solution
    white_pieces = ["White Pawn", "White Bishop"]
    black_pieces = []

    # Combine the lists and format the output string as per the requirements
    final_piece_list = white_pieces + black_pieces
    output_string = ", ".join(final_piece_list)

    print(output_string)

solve_chess_puzzle()