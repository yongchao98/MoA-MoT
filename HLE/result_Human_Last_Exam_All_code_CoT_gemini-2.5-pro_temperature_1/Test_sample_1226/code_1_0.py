def solve_chess_maximization():
    """
    Calculates the maximum material points for a mate in exactly 6 moves.

    The solution is based on a constructed chess position:
    - Black King (BK) is on h8.
    - White King (WK) is on a1.
    - A cage is formed by a White Queen on e6 (guards g8) and a White Pawn on g6 (guards h7).
    - The mate is delivered by another pawn, starting on h5, promoting on h8.
    - This requires a Queen on g1 to clear a path first.

    The mating sequence (M#6):
    1. Q(g1)g2
    2. Qg3
    3. Qg4  (The queen is now clear of the h-pawn's path)
    4. h6
    5. h7
    6. h8=Q#

    The script calculates the total material points by:
    1. Defining the points of the pieces in the mating 'engine'.
    2. Calculating how many squares are available to be filled by extra Queens.
       - A square is available if it's not used by the engine/path
         and a Queen placed there would not check the BK on h8.
    3. Summing the points of the engine and the extra Queens.
    """

    # --- 1. Define the 'engine' for the Mate in 6 ---

    # Value of the pieces in the M#6 engine
    engine_queen_1_val = 9  # Queen on e6
    engine_queen_2_val = 9  # Queen on g1
    engine_pawn_1_val = 1   # Pawn on g6
    engine_pawn_2_val = 1   # Pawn on h5
    engine_pieces_value = engine_queen_1_val + engine_queen_2_val + engine_pawn_1_val + engine_pawn_2_val

    # Squares occupied by the engine pieces (including kings)
    # Kings do not contribute to the material score.
    engine_squares = {'a1', 'h8', 'e6', 'g6', 'h5', 'g1'}

    # Path squares for the mating plan that must be empty at the start
    path_squares = {'g2', 'g3', 'g4', 'h6', 'h7'}
    
    committed_squares = engine_squares.union(path_squares)

    # --- 2. Calculate available squares for extra Queens ---
    
    all_squares = {f'{col}{row}' for col in 'abcdefgh' for row in '12345678'}
    
    # Squares from which a Queen would check the BK on h8
    forbidden_squares = {
        # Rank 8
        'a8', 'b8', 'c8', 'd8', 'e8', 'f8', 'g8',
        # File h
        'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'h7',
        # Diagonal a1-h8
        'a1', 'b2', 'c3', 'd4', 'e5', 'f6', 'g7'
    }

    # The set of all squares that cannot hold an extra queen
    unavailable_for_queens = committed_squares.union(forbidden_squares)
    
    num_total_squares = 64
    num_unavailable_squares = len(unavailable_for_queens)
    
    num_added_queens = num_total_squares - num_unavailable_squares
    queen_value = 9
    added_queens_value = num_added_queens * queen_value

    # --- 3. Calculate and print the final score ---

    total_value = engine_pieces_value + added_queens_value

    print("The greatest number of points is calculated as follows:")
    print("Value of the core pieces creating the mate-in-6 setup + Value of extra Queens placed on the board.\n")
    print("Value of engine pieces:")
    print(f"  Queen 1 (e6): {engine_queen_1_val} points")
    print(f"  Queen 2 (g1): {engine_queen_2_val} points")
    print(f"  Pawn 1 (g6): {engine_pawn_1_val} point")
    print(f"  Pawn 2 (h5): {engine_pawn_2_val} point")
    print(f"Total for engine pieces = {engine_pieces_value} points.\n")

    print("Number of extra Queens:")
    print(f"  There are {num_total_squares} squares on the board.")
    print(f"  {len(unavailable_for_queens)} squares are unavailable (used by the engine, its path, or are checking squares).")
    print(f"  Number of available squares for Queens = {num_total_squares} - {len(unavailable_for_queens)} = {num_added_queens}.\n")

    print("Final point calculation:")
    # The user requested to see each number in the final equation.
    print(f"{engine_pieces_value} + {num_added_queens} * {queen_value} = {engine_pieces_value} + {added_queens_value} = {total_value}")
    
solve_chess_maximization()
>>> 353