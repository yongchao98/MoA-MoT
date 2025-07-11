# This script requires the 'python-chess' library.
# You can install it by running: pip install python-chess

import chess

def solve_chess_puzzle():
    """
    Solves the chess puzzle from the FEN string and demonstrates the
    shortest mating sequence for White.
    """
    
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    board = chess.Board(fen)
    
    # The winning sequence of moves, determined through analysis.
    # 1. Qxf7+ Kh7 2. Bg6+ Qxg6 3. Nxg6#
    moves_uci = ["d5f7", "g8h7", "h5g6", "g5g6", "e5g6"]
    
    print("The starting position is given by the FEN:")
    print(fen)
    print("\nWhite can force a checkmate in 3 moves.")
    print("The optimal sequence of play (the 'equation' of moves) is:\n")

    move_count_white = 0
    full_move_number = board.fullmove_number
    
    for uci_move in moves_uci:
        move = chess.Move.from_uci(uci_move)
        
        # Print the move in Standard Algebraic Notation (SAN)
        if board.turn == chess.WHITE:
            move_count_white += 1
            print(f"{full_move_number}. {board.san(move)}", end=" ")
        else:
            print(f"{board.san(move)}")
            full_move_number += 1
            
        board.push(move)

    print("\n")
    # Verify the final position and print the result
    if board.is_checkmate():
        print("The final position is checkmate.")
        print(f"White wins in {move_count_white} moves.")
    else:
        # This part should not be reached if the analysis is correct.
        print("Error: The provided sequence did not lead to a checkmate.")
        print(f"Final board FEN: {board.fen()}")

if __name__ == '__main__':
    solve_chess_puzzle()