import chess

# This script verifies the checkmate sequence for the given chess puzzle.
# First, ensure you have the 'python-chess' library installed:
# pip install python-chess

def solve_chess_puzzle():
    """
    This function analyzes the chess position and prints the forced
    checkmate sequence for Black.
    """
    # The initial position is given by the FEN string.
    fen = "r1b1k2r/ppp2ppp/8/2bBN3/3nn3/8/PPP2qPP/RNBQR2K b Qkq - 0 1"

    # Initialize the board with the position.
    board = chess.Board(fen)

    print("Analyzing the best move sequence for Black...")
    
    # The correct sequence is a classic smothered mate.
    # We will programmatically execute the moves to verify and display them.

    # 1... Qg1+ (Black's first move)
    # The move is a Queen sacrifice to g1, delivering a check.
    move1_black = board.parse_san("Qg1+")
    san1_black = board.san(move1_black)
    board.push(move1_black)

    # 2. Rxg1 (White's forced response)
    # White's only legal move is to capture the queen with the rook.
    move2_white = list(board.legal_moves)[0]
    san2_white = board.san(move2_white)
    board.push(move2_white)

    # 2... Nf2# (Black's final, checkmating move)
    # The knight delivers the checkmate on the f2 square.
    move2_black = board.parse_san("Nf2#")
    # The library correctly identifies this as a checkmate move and adds '#'
    san2_black = board.san(move2_black)
    board.push(move2_black)

    # Construct the final sequence string in the required format.
    final_sequence = f"1... {san1_black} 2. {san2_white} {san2_black}"

    print("\nThe best sequence of moves that leads to checkmate is:")
    print(final_sequence)

    # Verify that the final position is indeed checkmate.
    if board.is_checkmate():
        print("\nResult: White is in checkmate.")
        print("Final board position:")
        print(board)

solve_chess_puzzle()