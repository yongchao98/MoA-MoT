import chess
import io

def solve_chess_puzzle():
    """
    This function sets up a specific chess position and verifies the
    correct mating sequence. It then prints each move of the sequence.
    
    To run this code, you need to install the python-chess library:
    pip install python-chess
    """
    # The FEN string for the chess position
    fen = "r1b1k2r/ppp2ppp/8/2bBN3/3nn3/8/PPP2qPP/RNBQR2K b Qkq - 0 1"
    
    # Set up the board with the given position
    board = chess.Board(fen)
    
    # The best mating sequence (Smothered Mate)
    # 1... Qg1+
    # 2. Rxg1 Nf2#
    moves = ["Qg1+", "Rxg1", "Nf2#"]
    
    print("The best sequence of moves leading to checkmate is:")
    
    # 1. Black's move
    move_san = moves[0]
    board.push_san(move_san)
    print(f"1... {move_san}")
    
    # 2. White's move
    # White only has one legal move, which is forced
    move_san = moves[1]
    board.push_san(move_san)
    print(f"2. {move_san}")

    # 2. Black's move
    move_san = moves[2]
    # The # symbol for mate is optional when pushing a move
    board.push_san(move_san.replace('#', ''))
    print(f"2... {move_san}")
    
    # Verify that the final position is indeed checkmate
    if board.is_checkmate():
        print("\nResult: White is in checkmate.")
    else:
        print("\nThis sequence did not result in a checkmate.")

solve_chess_puzzle()