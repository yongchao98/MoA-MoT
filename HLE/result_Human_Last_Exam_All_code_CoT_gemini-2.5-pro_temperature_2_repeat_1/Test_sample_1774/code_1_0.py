def solve_chess_puzzle():
    """
    Analyzes the given chess position and prints the mating sequence.
    """
    moves_to_mate = 2

    print("This is a classic mate-in-2 puzzle based on a queen sacrifice and a subsequent double check.")
    print("-" * 40)
    print(f"The solution is a checkmate in {moves_to_mate} moves.")
    print("Here is the winning sequence:")
    print("-" * 40)

    # Move 1
    white_move_1 = "Qxh7+"
    black_move_1 = "Kxh7"
    print(f"Move 1 (White): {white_move_1}")
    print("White sacrifices the Queen by capturing the pawn on h7, delivering a check.")
    print(f"Move 1 (Black): {black_move_1}")
    print("Black has no other choice but to capture the Queen with the King.")
    print("")

    # Move 2
    white_move_2 = "Ng6++"
    print(f"Move 2 (White): {white_move_2}")
    print("White moves the Knight from e5 to g6, delivering a devastating double check from both the Knight and the Bishop on d3.")
    print("\nAnalysis of the final position:")
    print("- The Black King on h7 is in check by the Knight on g6 and the Bishop on d3.")
    print("- Because it is a double check, the King must move.")
    print("- The escape squares g8 and h8 are both controlled by the Knight on g6.")
    print("- The King cannot move to any other square. This is checkmate.")

    print("-" * 40)
    print(f"Total moves for White to mate: {moves_to_mate}")


if __name__ == '__main__':
    solve_chess_puzzle()
