# First, you may need to install the required library. You can do this by running:
# pip install python-chess

import chess

def solve_chess_puzzle():
    """
    This function sets up the specified chess position and demonstrates
    the forced checkmate sequence for White.
    """
    # The board position is represented using Forsyth-Edwards Notation (FEN).
    # 'w' indicates that it is White's turn to move.
    fen = "rn4rk/pbppq1pp/1p2pb2/4N3/3PN3/3B4/PPP2PPP/R3K2R w KQ - 0 1"
    board = chess.Board(fen)

    print("Initial Position (White to move):")
    print(board)
    print("-" * 30)

    # This puzzle demonstrates a forced mate in 3 moves.
    # We will walk through one of the main variations.
    # The sequence is: 1. Nxf7+ Kxf7 2. Ng5+ Kg6 3. Qh7#

    # Move 1 for White: Nxf7+
    # The knight on e5 captures on f7, delivering a check.
    move1_w = "Nxf7+"
    board.push_san(move1_w)

    # Move 1 for Black: Kxf7
    # Black is forced to respond to the check. Taking the knight is one of two legal moves.
    move1_b = "Kxf7"
    board.push_san(move1_b)
    
    # The '1' represents the first full move in the sequence.
    print(f"1. {move1_w} {move1_b}")
    print(board)
    print("-" * 30)

    # Move 2 for White: Ng5+
    # The other knight moves to g5, delivering another check.
    move2_w = "Ng5+"
    board.push_san(move2_w)

    # Move 2 for Black: Kg6
    # The black king must move again. Kg6 is one of its options.
    move2_b = "Kg6"
    board.push_san(move2_b)
    
    # The '2' represents the second full move.
    print(f"2. {move2_w} {move2_b}")
    print(board)
    print("-" * 30)

    # Move 3 for White: Qh7#
    # The white queen delivers the final blow, resulting in checkmate.
    move3_w = "Qh7#"
    board.push_san(move3_w)
    
    # The '3' represents the final move leading to mate.
    print(f"3. {move3_w}")
    print(board)
    print("-" * 30)

    # Verify that the position is checkmate and print the final answer.
    if board.is_checkmate():
        print("Result: Confirmed Checkmate.")
        print("The number of moves for White to mate the Black king is 3.")
    else:
        print("Result: The sequence did not lead to a checkmate.")

if __name__ == '__main__':
    solve_chess_puzzle()
