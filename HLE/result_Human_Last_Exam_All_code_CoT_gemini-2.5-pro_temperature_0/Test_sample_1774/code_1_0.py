def solve_chess_puzzle():
    """
    This function prints the step-by-step solution to the described chess puzzle.
    """
    # The number of moves for White to deliver checkmate.
    moves_to_mate = 2

    # The move sequence for the checkmate.
    # Each element in the list represents one full move (White and Black).
    # The final move is just White's move.
    sequence = [
        {"number": 1, "white_move": "Qxh7+", "black_move": "Kxh7"},
        {"number": 2, "white_move": "Ng6#", "black_move": ""}
    ]

    print(f"White can checkmate the Black king in {moves_to_mate} moves.")
    print("The winning sequence is a queen sacrifice followed by a double check:")
    print("-" * 30)

    # Print each move in the sequence
    for move in sequence:
        # The final move doesn't have a response from Black
        if move["black_move"]:
            print(f"Move {move['number']}: {move['white_move']} {move['black_move']}")
        else:
            print(f"Move {move['number']}: {move['white_move']}")

    print("-" * 30)
    print("Explanation:")
    print("1. Qxh7+: White sacrifices the queen. This is a check, and Black's only legal move is to capture the queen with the king.")
    print("2. Ng6#: The knight moves to G6, delivering a double check from both the knight and the bishop on D3. The Black king has no squares to move to, resulting in checkmate.")

if __name__ == "__main__":
    solve_chess_puzzle()
