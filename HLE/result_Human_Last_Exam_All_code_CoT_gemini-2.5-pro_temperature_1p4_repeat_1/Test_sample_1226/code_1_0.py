def calculate_max_material():
    """
    Calculates the maximum material value for a mate in exactly 6 moves.

    The calculation is based on a specific chess position:
    - A Rook at a2 travels for 5 moves to set up a mate on the 6th move.
    - The mate is Raxg7#, capturing a friendly Queen on g7.
    - The Black King is caged by pairs of Queens that block each other.
    - All remaining squares are filled with filler Queens.
    """
    
    # Standard piece values
    queen_val = 9
    rook_val = 5
    
    # --- Core Pieces for the Mate and Cage ---
    
    # 1. The attacking rook
    attacking_rook = rook_val
    print(f"Attacking Rook: {attacking_rook} points")
    
    # 2. The piece to be captured on g7 to enable the mate. Make it a Queen for max value.
    sacrificial_queen = queen_val
    print(f"Sacrificial Queen on g7: {sacrificial_queen} points")
    
    # 3. Three pairs of Queens used to form the cage. Each pair consists of one
    #    caging Queen and one blocking Queen. (3 * 2 = 6 Queens)
    caging_queens_count = 6
    caging_queens_value = caging_queens_count * queen_val
    print(f"Caging and Blocking Queens (6 pieces): {caging_queens_value} points")
    
    core_material_value = attacking_rook + sacrificial_queen + caging_queens_value
    print("---------------------------------")
    print(f"Subtotal for Core Pieces: {core_material_value} points")
    print("---------------------------------")

    # --- Filler Pieces ---
    
    # Calculate the number of squares used by the setup
    total_squares = 64
    
    # Kings
    king_squares = 2 # WK, BK
    
    # Core pieces from above
    core_piece_squares = 1 + 1 + caging_queens_count # Rook, Sacrificial Queen, 6 caging Queens
    
    # Squares for the Rook's path (a3, a4, a5, a6, a7)
    rook_path_squares = 5
    
    # Square for the Black King's forced move (g8)
    black_king_move_square = 1
    
    used_squares = king_squares + core_piece_squares + rook_path_squares + black_king_move_square
    
    filler_squares = total_squares - used_squares
    filler_queens_value = filler_squares * queen_val
    print(f"Filler Queens on remaining squares ({filler_squares} pieces): {filler_queens_value} points")

    # --- Total Value ---
    total_value = core_material_value + filler_queens_value
    print("=================================")
    print(f"Total material value = {core_material_value} (Core) + {filler_queens_value} (Filler)")
    print(f"Final Answer: {total_value}")
    
    return total_value

if __name__ == '__main__':
    calculate_max_material()