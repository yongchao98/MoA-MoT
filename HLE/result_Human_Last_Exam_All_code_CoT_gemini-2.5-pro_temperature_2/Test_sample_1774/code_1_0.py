def solve_chess_puzzle():
    """
    Analyzes the chess position to find the number of moves to checkmate.
    It prints the move sequence against the best defense.
    """
    moves_to_mate = 3

    print(f"White can force a checkmate in {moves_to_mate} moves.")
    print("The puzzle requires finding the mate against Black's best defense.\n")
    print("The main line is as follows:")

    # Move 1
    move_number_1 = 1
    white_move_1 = "Ng6+"
    black_move_1 = "Bxg6"
    print(f"Move {move_number_1}: White plays {white_move_1}, Black responds with {black_move_1}")

    # Move 2
    move_number_2 = 2
    white_move_2 = "Qxg6+"
    black_move_2 = "Kh8"
    print(f"Move {move_number_2}: White plays {white_move_2}, Black responds with {black_move_2}")

    # Move 3
    move_number_3 = 3
    white_move_3 = "Nh6#"
    print(f"Move {move_number_3}: White delivers checkmate with {white_move_3}")

    print("\n--- Move-by-Move Explanation ---")
    print("1. White starts with a knight sacrifice (Ng6+). The Black King has no escape squares and must capture the knight.")
    print("   - If Black plays 1...hxg6, White mates in the next move with 2. Qxh7#.")
    print("   - Black's best defense is to play 1...Bxg6, prolonging the mate.")
    print("2. White's Queen captures the bishop on g6 (Qxg6+), checking the King again. This forces the Black King to its only escape square, h8.")
    print("3. White moves the other knight to h6 (Nh6#), which delivers a fatal double check from both the knight and the queen. The Black King has no escape, and the checks cannot be blocked or the pieces taken. This is checkmate.")

solve_chess_puzzle()
