def solve_diagonal_corridor_mate():
    """
    This function determines the minimum pieces required for a Diagonal Corridor Mate
    under the specified constraints and prints the result.
    """

    # Piece values are used for sorting and for calculating material advantage.
    piece_values = {
        'Pawn': 1,
        'Knight': 3,
        'Bishop': 3,
        'Rook': 5,
        'Queen': 9
    }

    # Based on the step-by-step logical deduction:
    # 1. White needs a Bishop (minimal value) to deliver the diagonal check.
    # 2. Black needs pieces on g8 and h7. To satisfy the material advantage
    #    constraint (Black > White), Black's pieces must have a total value
    #    greater than the White Bishop's value (3).
    # 3. The lowest value combination for Black is a Pawn (1) and a Knight (3),
    #    totaling 4. A Pawn and Bishop is an equivalent alternative. We choose one.
    additional_pieces = [
        {'color': 'White', 'type': 'Bishop'},
        {'color': 'Black', 'type': 'Pawn'},
        {'color': 'Black', 'type': 'Knight'}
    ]

    # Separate the pieces by color
    white_pieces = [p for p in additional_pieces if p['color'] == 'White']
    black_pieces = [p for p in additional_pieces if p['color'] == 'Black']

    # Sort each list by piece value in ascending order
    white_pieces.sort(key=lambda p: piece_values[p['type']])
    black_pieces.sort(key=lambda p: piece_values[p['type']])

    # Format the names for the final output string
    white_piece_names = [f"{p['color']} {p['type']}" for p in white_pieces]
    black_piece_names = [f"{p['color']} {p['type']}" for p in black_pieces]

    # Combine the lists and join into a comma-separated string
    final_list = white_piece_names + black_piece_names
    result_string = ", ".join(final_list)

    print(result_string)

solve_diagonal_corridor_mate()