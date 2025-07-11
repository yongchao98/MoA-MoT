import chess

def solve_chess_puzzle():
    """
    This function solves the given chess puzzle by executing the correct
    sequence of moves and printing the result.
    """
    # The initial position in Forsyth-Edwards Notation (FEN)
    fen = "r1b1k2r/ppp2ppp/8/2bBN3/3nn3/8/PPP2qPP/RNBQR2K b Qkq - 0 1"
    board = chess.Board(fen)

    print("Initial Position (FEN):")
    print(board.fen())
    print(board)
    print("\nExecuting the winning sequence for Black...\n")

    # The correct sequence of moves from option C
    moves = ["Qg1+", "Nf2#"]
    
    # 1. Black's move
    move_uci = board.push_san(moves[0])
    print(f"1. ... {moves[0]}")
    print(board)
    print(f"Move in UCI notation: {move_uci}\n")

    # 2. White's move
    # The library will automatically find the only legal move for white
    white_move = board.find_move(board.turn, chess.KING)
    # The python-chess library will automatically determine the only legal move for White
    white_move_san = board.san(board.legal_moves.pop())
    board.push(white_move_san)
    print(f"2. {white_move_san}", end=" ")
    
    # 2. Black's second move
    move_uci = board.push_san(moves[1])
    print(f"{moves[1]}")
    print(board)
    print(f"Move in UCI notation: {move_uci}\n")

    # Verify if it's checkmate
    if board.is_checkmate():
        print("Final Result: CHECKMATE! Black wins.")
    else:
        print("Final Result: Not a checkmate.")

if __name__ == '__main__':
    # For this script to run, you need to install the 'python-chess' library:
    # pip install python-chess
    solve_chess_puzzle()
