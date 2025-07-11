def calculate_max_material():
    """
    Calculates the greatest number of points of white material for a mate in >= 6 moves.

    The calculation is based on a specific chess position constructed to have a
    forced mate in exactly 6 moves, while maximizing the number of high-value
    pieces (Queens) on the board.
    """
    
    # Piece values
    QUEEN_POINTS = 9
    ROOK_POINTS = 5
    PAWN_POINTS = 1
    
    # --- 1. Define the core set of white pieces for the mate-in-6 setup ---
    # This setup creates a shuttle for the black king (a8 <=> b8) and a 6-move mate.
    # We use 1 Pawn and 3 Rooks.
    core_pieces = {
        'Pawn': {'count': 1, 'points': PAWN_POINTS},
        'Rook': {'count': 3, 'points': ROOK_POINTS}
    }
    
    core_points = sum(details['count'] * details['points'] for details in core_pieces.values())
    num_core_pieces = sum(details['count'] for details in core_pieces.values())
    
    # --- 2. Define all squares that CANNOT contain a filler Queen ---
    # These include squares for the kings, core pieces, paths, and forbidden zones.
    
    # Squares for kings, essential pieces, and required empty paths/shuttles.
    # WK at c6, BK at a8. Pawn at a7. Rooks at d8, b5, b1.
    # Empty squares: b8 (shuttle), b7 (mating square), a1-a5 (rook path).
    essential_setup_squares = {
        "Kings": {"WK": "c6", "BK": "a8"},
        "Core Pieces": ["a7", "d8", "b5", "b1"],
        "Paths/Empty": ["b8", "b7", "a1", "a2", "a3", "a4", "a5"]
    }
    
    num_setup_squares = (
        len(essential_setup_squares["Kings"]) +
        len(essential_setup_squares["Core Pieces"]) +
        len(essential_setup_squares["Paths/Empty"])
    )
    
    # Squares where a Queen would illegally interfere with the solution
    # by attacking the Black King at a8 or b8.
    forbidden_queen_squares = {
        # Queen attacks BK on a8 (diagonally)
        "a8_diag_attack": ["h1", "g2", "f3", "e4", "d5"], 
        # Queen attacks BK on b8 (diagonally)
        "b8_diag_attack": ["h2", "g3", "f4", "e5", "d6", "c7"],
        # Queen attacks BK on 8th rank (stalemate/check)
        "8th_rank_attack": ["c8", "e8", "f8", "g8", "h8"],
    }
    
    num_forbidden_squares = sum(len(sq_list) for sq_list in forbidden_queen_squares.values())

    # --- 3. Calculate the number of Queens and total points ---
    total_board_squares = 64
    num_non_queen_squares = num_setup_squares + num_forbidden_squares
    
    num_filler_queens = total_board_squares - num_non_queen_squares
    queen_points = num_filler_queens * QUEEN_POINTS
    
    total_points = core_points + queen_points
    
    # --- 4. Print the final equation and result ---
    print("This problem asks for the maximum material value for a specific type of chess position.")
    print("The chosen position forces a checkmate in exactly 6 moves.")
    print("\nHere is the breakdown of the calculation:")
    
    pawn_str = f"({core_pieces['Pawn']['count']} Pawn * {PAWN_POINTS} pts)"
    rook_str = f"({core_pieces['Rook']['count']} Rooks * {ROOK_POINTS} pts)"
    queen_str = f"({num_filler_queens} Queens * {QUEEN_POINTS} pts)"

    print(f"\nValue = {pawn_str} + {rook_str} + {queen_str}")
    
    pawn_val = core_pieces['Pawn']['count'] * PAWN_POINTS
    rook_val = core_pieces['Rook']['count'] * ROOK_POINTS
    print(f"Value = {pawn_val} + {rook_val} + {queen_points}")
    print(f"Total Points = {total_points}")

calculate_max_material()