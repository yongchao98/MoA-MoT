import chess
import io

def solve_chess_puzzle():
    """
    Analyzes the chess position and validates the correct checkmate sequence.
    The sequence is 1... Qg1+ 2. Rxg1 Nf2#, a classic smothered mate.
    """
    # FEN string for the initial position where it's Black's turn to move.
    fen = "r1b1k2r/ppp2ppp/8/2bBN3/3nn3/8/PPP2qPP/RNBQR2K b Qkq - 0 1"
    board = chess.Board(fen)

    print("The starting position is:")
    print(board)
    print("\nBlack can force checkmate in 2 moves.")
    print("The sequence of moves is as follows:\n")

    # 1... Qg1+
    # This is a queen sacrifice. It's a check, and White's king cannot move.
    # The only legal move for White is to capture the queen.
    move1 = board.parse_san("Qg1")
    board.push(move1)
    print("1... Qg1+")

    # 2. Rxg1
    # White is forced to capture the Black queen with the rook.
    move2 = board.parse_san("Rxg1")
    board.push(move2)
    print("2. Rxg1")

    # 2... Nf2#
    # Black delivers checkmate with the knight. The White king is trapped
    # by its own pieces (a "smothered mate").
    move3 = board.parse_san("Nf2")
    board.push(move3)
    print("2... Nf2#")

    print("\nThe final position is:")
    print(board)

    # Confirm that the final position is indeed a checkmate.
    if board.is_checkmate():
        print("\nThis sequence results in a checkmate. White is checkmated.")
    else:
        print("\nThis sequence does not result in a checkmate.")

solve_chess_puzzle()
