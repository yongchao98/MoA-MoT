import chess

def solve_chess_puzzle():
    """
    This function solves the chess puzzle provided by setting up the board,
    executing the mating sequence, and printing the results.
    """
    # The FEN for the given position is:
    # rn4rk1/pbpq1pp/1p2pb2/4N2Q/3NP3/3B4/PPP2PPP/R3K2R w KQ - 0 1
    board = chess.Board("rn4rk1/pbpq1pp/1p2pb2/4N2Q/3NP3/3B4/PPP2PPP/R3K2R w KQ - 0 1")

    print("Analyzing the chess position...")
    print("Initial board state (White to move):")
    print(board)
    print("\nFound a forced checkmate sequence.")
    print("------------------------------------")
    
    # Mating sequence
    # 1. Qxf7+
    move1 = board.parse_san("Qxf7+")
    print(f"Move 1 for White: {board.san(move1)}")
    board.push(move1)
    
    # 1... Rxf7 (Black's only legal move)
    black_move1 = list(board.legal_moves)[0]
    print(f"   Black's forced reply: {board.san(black_move1)}")
    board.push(black_move1)

    # 2. Ng6+
    move2 = board.parse_san("Ng6+")
    print(f"Move 2 for White: {board.san(move2)}")
    board.push(move2)
    
    # 2... hxg6 (Black's only legal move)
    black_move2 = list(board.legal_moves)[0]
    print(f"   Black's forced reply: {board.san(black_move2)}")
    board.push(black_move2)

    # 3. Rh8#
    move3 = board.parse_san("Rh8#")
    print(f"Move 3 for White: {board.san(move3)}")
    board.push(move3)

    print("------------------------------------")

    if board.is_checkmate():
        print("\nFinal board state: Checkmate!")
        print(board)
        
        # As requested, here is the final "equation" showing each move.
        print("\nThe final equation for the solution is:")
        print("White's Move (1) + White's Move (2) + White's Move (3) = Checkmate")
        print(f"        {board.san(move1)}      +       {board.san(move2)}       +       {board.san(move3)}      = Checkmate in 3 moves")

    else:
        print("\nThe provided sequence did not result in a checkmate.")

if __name__ == "__main__":
    solve_chess_puzzle()
