# To run this code, you may need to install the python-chess library:
# pip install python-chess

import chess

def demonstrate_winning_move():
    """
    Sets up the chess position from the FEN and demonstrates the winning line for White.
    """
    # The FEN string for the given chess position.
    fen = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"
    board = chess.Board(fen)

    print("--- Initial Position ---")
    print("White has a winning advantage due to the a7-pawn.")
    print(board)
    print(f"FEN: {board.fen()}\n")

    # The winning move is Nc5. The best reply for black is Nxc5.
    # The sequence demonstrates how white achieves a winning endgame.
    # This sequence is our "final equation" as requested, with numbers 1, 5, 5, 2, 5 from the moves.
    print("--- The Winning Equation ---")
    print("The best move is 1. Nc5. If Black plays the best response, the sequence is:")
    print("1. Nc5 Nxc5 2. bxc5\n")

    # --- Move 1: White plays Nc5 ---
    move_w1 = board.parse_san("Nc5")
    board.push(move_w1)
    print("--- Position after 1. Nc5 ---")
    print("White attacks the knight's defender and forces a reaction.")
    print(board)
    print(f"FEN: {board.fen()}\n")

    # --- Move 1: Black plays Nxc5 ---
    move_b1 = board.parse_san("Nxc5")
    board.push(move_b1)
    print("--- Position after 1... Nxc5 ---")
    print("Black is forced to trade knights.")
    print(board)
    print(f"FEN: {board.fen()}\n")

    # --- Move 2: White plays bxc5 ---
    move_w2 = board.parse_san("bxc5")
    board.push(move_w2)
    print("--- Position after 2. bxc5 ---")
    print("White creates a second passed pawn. Black cannot stop both pawns from promoting.")
    print(board)
    print(f"FEN: {board.fen()}\n")

if __name__ == '__main__':
    demonstrate_winning_move()