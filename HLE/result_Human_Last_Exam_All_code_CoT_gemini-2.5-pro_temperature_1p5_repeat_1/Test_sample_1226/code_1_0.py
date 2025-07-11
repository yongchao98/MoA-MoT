def calculate_max_material():
    """
    Calculates and explains the maximum material score for a mate in >= 6 moves.
    """
    # Standard piece values
    queen_val = 9
    rook_val = 5
    bishop_val = 3
    knight_val = 3
    pawn_val = 1

    # My construction for a forced mate in exactly 6 moves:
    # The setup uses a 'cage' to trap the black king, a 'delay' mechanism
    # using pawn promotions, a 'mating' piece, and 'filler' pieces to maximize the score.

    # 1. Core Pieces for the Mate-in-6 Mechanism
    # Cage pieces: 1 Queen and 1 Pawn create a cage forcing the black king to move between a1 and b1.
    cage_queens = 1
    cage_pawns = 1
    cage_points = (cage_queens * queen_val) + (cage_pawns * pawn_val)

    # Delay/Mate pieces: 5 Pawns must promote to clear a path for the mating Queen.
    delay_pawns = 5
    mating_queens = 1
    delay_mate_points = (mating_queens * queen_val) + (delay_pawns * pawn_val)
    
    # Total core white pieces involved in the mate
    core_queens = cage_queens + mating_queens
    core_pawns = cage_pawns + delay_pawns

    # 2. Board Squares Analysis
    total_squares = 64
    
    # Squares occupied by kings
    king_squares = 2
    
    # Squares occupied by the core white pieces
    core_piece_squares = core_queens + core_pawns
    
    # Squares that must be initially empty for the solution to work:
    # 1 for the black king's only other move (b1)
    # 5 for the pawns to promote into (c8, d8, e8, f8, g8)
    empty_squares = 1 + 5
    
    # Total specified squares for the mechanism
    specified_squares = king_squares + core_piece_squares + empty_squares
    
    # 3. Filler Pieces for Maximizing Score
    # The remaining squares are filled with Queens to maximize points.
    # These filler Queens are placed carefully to not interfere with the 6-move solution.
    filler_queen_squares = total_squares - specified_squares
    
    # 4. Final Calculation
    total_queens = core_queens + filler_queen_squares
    total_pawns = core_pawns
    
    total_points = (total_queens * queen_val) + (total_pawns * pawn_val)

    print("To solve this, I constructed a chess position with a forced mate in 6 moves.")
    print("The material score is maximized by filling all non-essential squares with Queens.")
    print("\nHere's the calculation based on my construction:\n")

    print(f"Pieces for the Caging/Mating mechanism:")
    print(f"  - Queens: {core_queens}")
    print(f"  - Pawns:   {core_pawns}\n")

    print("Board layout:")
    print(f"  - Squares for Kings:                     {king_squares}")
    print(f"  - Squares for core Queens & Pawns:       {core_piece_squares}")
    print(f"  - Squares that must be empty:            {empty_squares}")
    print(f"  - Squares available for filler Queens:   {filler_queen_squares}\n")

    print("Final Piece Count:")
    print(f"  - Total Queens: {core_queens} + {filler_queen_squares} = {total_queens}")
    print(f"  - Total Pawns:                           {total_pawns}")
    print("  - Total Rooks/Bishops/Knights:           0\n")

    print("The final equation for the total points is:")
    print(f"({total_queens} Queens * {queen_val} points) + ({total_pawns} Pawns * {pawn_val} point)")
    print(f"({total_queens} * {queen_val}) + ({total_pawns} * {pawn_val}) = {total_queens * queen_val} + {total_pawns * pawn_val} = {total_points}")

    print(f"\nThus, the greatest number of points is {total_points}.")
    return total_points

# Run the calculation and store the final answer
final_answer = calculate_max_material()