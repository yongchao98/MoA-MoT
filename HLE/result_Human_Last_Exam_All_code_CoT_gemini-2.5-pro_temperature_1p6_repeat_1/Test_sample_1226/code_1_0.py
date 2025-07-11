def solve_chess_puzzle():
    """
    Calculates the maximum material points for a mate in exactly 6 moves.

    The solution is based on a constructed chess position where:
    1. The black king is trapped and can only shuttle between two squares.
    2. A line of 5 white "blocking" queens must be moved one by one (5 moves) to clear a path for a 6th "mating" queen.
    3. The 6th move is the checkmate.
    4. The rest of the board is filled with queens to maximize points.
    """
    
    # Standard material point values
    piece_values = {
        'Queen': 9,
        'Rook': 5,
        'Bishop': 3,
        'Knight': 3,
        'Pawn': 1
    }

    # Piece counts based on the constructed mate-in-6 position
    
    # Containment pieces to trap the black king
    # 3 Pawns (on h6, f6, g5) and 1 Queen (on f1)
    containment_queens = 1
    containment_pawns = 3
    
    # Mating mechanism pieces
    # 1 Mating Queen and 5 Blocking Queens on the 7th rank
    mating_and_blocking_queens = 6
    
    # Filler pieces
    # The board has 64 squares. 2 are for kings.
    # 4 squares are for containment pieces.
    # 6 squares are for mating/blocking pieces.
    # Remaining squares = 64 - 2 - 4 - 6 = 52.
    # We fill these with the highest value piece, the Queen.
    filler_queens = 52
    
    # Total piece counts
    total_queens = containment_queens + mating_and_blocking_queens + filler_queens
    total_pawns = containment_pawns
    total_rooks = 0
    total_bishops = 0
    total_knights = 0
    
    # Calculate the total material points
    total_points = (total_queens * piece_values['Queen']) + \
                   (total_pawns * piece_values['Pawn'])

    print("To achieve a forced mate in exactly 6 moves with maximum material, we construct a specific board state.")
    print("The configuration uses the following white pieces:")
    print(f"- Queens: {total_queens}")
    print(f"- Pawns: {total_pawns}")
    print("\nThe calculation for the total material points is as follows:")
    
    # Output the final equation and result
    queen_val = piece_values['Queen']
    pawn_val = piece_values['Pawn']
    
    print(f"{total_queens} Queens * {queen_val} points + {total_pawns} Pawns * {pawn_val} point = {total_points} points")

solve_chess_puzzle()
<<<534>>>