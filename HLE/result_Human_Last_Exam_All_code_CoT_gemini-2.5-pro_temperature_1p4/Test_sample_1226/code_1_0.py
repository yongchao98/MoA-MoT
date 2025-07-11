def solve_chess_puzzle():
    """
    Calculates the maximum white material for a mate in exactly 6 moves.

    The solution is constructed based on a common chess problem theme where a mate is
    delayed by White's own pieces, and Black is given tempo moves by a pawn.

    The constructed position and mating plan are as follows:
    - Black King on a8.
    - A stalemating trap is formed by a White King on c8 and a White Queen on c7.
    - A Black Pawn on h7 provides 5 tempo moves for Black (...h6, ...h5, etc.).
    - White's mating plan is to play a Queen to b7. This Queen starts at b1.
    - The path for the mating Queen is blocked by 5 other White Queens on b2, b3, b4, b5, and b6.
    - White's first 5 moves are used to move these blocking Queens away. This requires
      at least one empty square to start the clearing sequence. For subsequent moves,
      the Queens can move into the squares vacated by the advancing Black Pawn.
    - The 6th move is the checkmate (Qb7#).
    - All other available squares are filled with White Queens to maximize the material score.
    """

    # --- Setup and Calculation ---

    QUEEN_VAL = 9
    total_squares = 64

    # Squares occupied by pieces other than White Queens, or are empty:
    # 1. White King (essential for the trap)
    # 2. Black King (the target)
    # 3. Black Pawn (provides tempo for Black's 5 moves)
    # 4. A single empty square (required for the first blocking queen to move to)
    num_wk = 1
    num_bk = 1
    num_bp = 1
    num_empty = 1

    non_queen_squares = num_wk + num_bk + num_bp + num_empty

    # The remaining squares can all be filled with White Queens.
    num_white_queens = total_squares - non_queen_squares

    # Calculate the total material value
    total_material_value = num_white_queens * QUEEN_VAL

    # --- Output the results ---
    print("Plan: Calculate max material for a mate in exactly 6 moves.")
    print("-" * 60)
    print("Configuration:")
    print(" - Black King at a8, stalemated by White King at c8 & White Queen at c7.")
    print(" - Black Pawn at h7 provides 5 tempo moves for Black.")
    print(" - White's mating move (from b1 to b7) is blocked by 5 Queens.")
    print(" - White takes 5 moves to clear the path, requiring one empty square.")
    print(" - The 6th move is checkmate.")
    print("-" * 60)
    print("Calculation of Total White Queens:")
    print(f"  Total squares on board: {total_squares}")
    print(f"  Squares used by non-queen pieces or empty:")
    print(f"    - White King: {num_wk}")
    print(f"    - Black King: {num_bk}")
    print(f"    - Black Pawn: {num_bp}")
    print(f"    - Empty Squares: {num_empty}")
    print(f"  Total unavailable squares: {non_queen_squares}")
    print(f"  Number of White Queens = {total_squares} - {non_queen_squares} = {num_white_queens}")
    print("-" * 60)
    print("Final Material Value Calculation:")
    print(f"  Total White Material = (Number of White Queens) * (Value of a Queen)")
    print(f"                       = {num_white_queens} * {QUEEN_VAL} = {total_material_value}")
    print("-" * 60)
    print(f"The greatest number of points of white material is: {total_material_value}")

if __name__ == "__main__":
    solve_chess_puzzle()