def solve_chess_puzzle():
    """
    Calculates the greatest number of points of white material for a mate-in-6 position.
    
    The position is constructed as follows:
    1. A core mate-in-6 mechanism on the g and h files:
       - White King on g6, White Rook on h1.
       - Black King on h8 (shuttling to g8).
       - White plays 1.Rh2 ... 6.Rh7#
    
    2. A "prison" of immobilized pieces on files a-f to maximize points:
       - A 6x8 grid of interlocking Queens and Pawns that have no legal moves.
       - Rows 1, 3, 5, 7 (files a-f) are Queens.
       - Rows 2, 4, 6, 8 (files a-f) are Pawns.
    """
    
    # Define piece values
    queen_value = 9
    rook_value = 5
    pawn_value = 1
    
    # Count the pieces in our construction
    # The 'prison' on files a-f
    num_queens_in_prison = 6 * 4  # 6 columns, 4 rows of Queens
    num_pawns_in_prison = 6 * 4   # 6 columns, 4 rows of Pawns
    
    # The piece for the mating attack
    num_rooks_for_mate = 1
    
    # Calculate total points
    queen_points = num_queens_in_prison * queen_value
    pawn_points = num_pawns_in_prison * pawn_value
    rook_points = num_rooks_for_mate * rook_value
    
    total_points = queen_points + pawn_points + rook_points
    
    # Output the explanation and calculation
    print("To achieve the maximum material for a mate-in-6, we construct a board with a specific mating pattern and a large block of immobilized high-value pieces.")
    print("The final calculation is based on the following white pieces:")
    print(f"- {num_queens_in_prison} Queens from the immobilized block.")
    print(f"- {num_pawns_in_prison} Pawns from the immobilized block.")
    print(f"- {num_rooks_for_mate} Rook to perform the checkmate.")
    print("\nFinal point calculation:")
    print(f"{num_queens_in_prison} Queens * {queen_value} points/Queen = {queen_points} points")
    print(f"{num_pawns_in_prison} Pawns * {pawn_value} point/Pawn = {pawn_points} points")
    print(f"{num_rooks_for_mate} Rook * {rook_value} points/Rook = {rook_points} points")
    print(f"\nTotal points = {queen_points} + {pawn_points} + {rook_points} = {total_points}")

solve_chess_puzzle()