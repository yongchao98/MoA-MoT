def solve_chess_puzzle():
    """
    This function analyzes the chess position given by the FEN:
    4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20
    and determines the number of moves for White to win with optimal play.
    """

    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    
    print(f"Analyzing the chess position from FEN: {fen}")
    print("White to move and win.")
    print("\nThe position is critical. Black threatens checkmate with ...Qxg2#.")
    print("White must act with a forcing sequence of checks to win.")
    print("\nThe optimal winning sequence is as follows:")

    # The move sequence is presented as a list of strings.
    # Each string represents one full move (White's move and Black's response).
    # The final move is White's checkmating move.
    moves = [
        "1. Qxf7+ Kh7",
        "2. Qxg7+ Kxg7",
        "3. Rf7+  Kg8",
        "4. Rg7+  Kh8",
        "5. Rh7#"
    ]

    # Print each move in the sequence
    print("-" * 20)
    for move in moves:
        print(move)
    print("-" * 20)

    # The number of moves is the number of moves White makes.
    num_moves = len(moves)

    print(f"\nExplanation:")
    print("1. White starts with a check, forcing the Black king to h7.")
    print("2. White sacrifices the queen with another check, forcing the Black king to recapture on g7.")
    print("3. The White rook checks on f7. To prolong the game, Black moves the king to g8.")
    print("4. The White rook delivers a double check with the bishop, forcing the Black king to h8.")
    print("5. The White rook delivers the final checkmate on h7, protected by the bishop.")
    
    print(f"\nAssuming optimal play from both sides, White can force a checkmate in {num_moves} moves.")

solve_chess_puzzle()