def solve_chess_puzzle():
    """
    This function explains the solution to the chess puzzle and calculates the moves to mate.
    """

    # The mating sequence involves two distinct moves from White.
    move_1 = 1  # Queen sacrifice
    move_2 = 1  # Knight and Bishop deliver the checkmate

    total_moves_to_mate = move_1 + move_2

    print("White has a forced checkmate in 2 moves.")
    print("-----------------------------------------")
    print("The winning sequence is:")
    print("1. Qxh7+  (White Queen captures the pawn on H7, checking the Black King)")
    print("   Black's only move is Kxh7 (King captures the Queen).")
    print("\n2. Ng6#   (White Knight moves to G6, delivering a decisive double check)")
    print("   This is checkmate, as the Black King has no escape squares.")
    print("-----------------------------------------")
    print("The total number of moves for White in this sequence is calculated as follows:")
    
    # Printing each number in the final equation as requested.
    print(f"{move_1} + {move_2} = {total_moves_to_mate}")

if __name__ == "__main__":
    solve_chess_puzzle()