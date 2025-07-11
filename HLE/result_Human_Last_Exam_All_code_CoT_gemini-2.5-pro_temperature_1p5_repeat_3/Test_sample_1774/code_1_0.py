def solve_chess_puzzle():
    """
    This function explains and prints the solution to the chess puzzle.
    The solution is a forced checkmate in 4 moves.
    """

    print("The starting position is complex, but White has a forced checkmate.")
    print("The key is to use the knights to break open Black's kingside defenses, allowing the queen to deliver the final blows.")
    print("Black will attempt to prolong the game, so we must find the line that checkmates against the best defense.")
    print("\nThe solution is a mate in 4 moves. Here is the sequence:\n")

    move_1_white = "1. N(e5)xg6+"
    move_1_black = "Kh8"
    explanation_1 = "White sacrifices the knight on g6. If Black plays hxg6, it's mate in 2 with 2. Qxh7#. If Black plays fxg6, it's mate in 3 (2. Nxf6+ Kh8 3. Qxh7#). Black's best defense is moving the King."

    move_2_white = "2. Qxh7+"
    move_2_black = "Kxh7"
    explanation_2 = "White sacrifices the Queen. The pawn on g7 cannot capture because it is pinned to the king by the White knight on g6. The King is forced to capture."

    move_3_white = "3. N(g6)f8++"
    move_3_black = "Kg8"
    explanation_3 = "This is a double check from the Knight on f8 and a discovered check from the Rook on h1. The Black King has only one legal move."

    move_4_white = "4. Rh8#"
    explanation_4 = "The Rook on h8 delivers checkmate. The Black King cannot escape to f7 (covered by the Knight on e4) or f8 (covered by the Knight on f8)."

    final_answer = 4

    print(f"{move_1_white:<15} {move_1_black:<15} | {explanation_1}")
    print(f"{move_2_white:<15} {move_2_black:<15} | {explanation_2}")
    print(f"{move_3_white:<15} {move_3_black:<15} | {explanation_3}")
    print(f"{move_4_white:<15} {'':<15} | {explanation_4}")

    print(f"\nThus, the number of moves to mate the black king is {final_answer}.")


solve_chess_puzzle()