# First, ensure you have the required library installed:
# pip install python-chess

import chess

def find_mate_in_two():
    """
    This function sets up the described chess position and executes
    the moves for a mate in 2, printing the solution.
    """
    # FEN representing the board state:
    # Black: R(a8),N(b8),R(f8),K(g8); p(a7,c7,d7,g7,h7,b6,e6); B(b7,f6); Q(e7)
    # White: N(e4,e5); Q(h5); p(d4); B(d3); p(a2,b2,c2,f2,g2,h2); R(a1,h1); K(e1)
    # White to move.
    fen = "rn2r1k/pbppq1pp/1p2pb2/4N2Q/3PN3/3B4/PPP2PPP/R3K2R w KQkq - 0 1"
    
    board = chess.Board(fen)

    print("--- Chess Puzzle: White to Mate ---")
    print("\nInitial Board Position:")
    print(board)
    
    # --- The Mating Sequence ---
    # 1. Qxh7+ : The Queen, supported by the Bishop on d3, captures the h7 pawn,
    #            checking the King on g8. Black has only one legal move.
    # 1... Kxh7: The King is forced to capture the Queen.
    # 2. Ng6#  : The Knight from e5 moves to g6, delivering a check. The King on h7
    #            cannot move to g8 or h8 (attacked by the Knight) and cannot
    #            capture the Knight because it is defended by the second Knight on e4.
    #            This is checkmate.

    # We will now execute and verify these moves.
    try:
        # 1. White's first move
        move1_w_san = "Qxh7+"
        board.push_san(move1_w_san)

        # 1. Black's forced response
        move1_b_san = "Kxh7"
        board.push_san(move1_b_san)
        
        # 2. White's second move (the checkmate)
        # We find the Standard Algebraic Notation for the move from e5 to g6.
        move2_w_obj = chess.Move.from_uci("e5g6")
        move2_w_san = board.san(move2_w_obj)
        board.push(move2_w_obj)

        if board.is_checkmate():
            print("\nSolution Found: Mate in 2")
            print("The final move results in a checkmate.")
            
            print("\nFinal Equation:")
            # Output each number in the final equation
            print(f"1. {move1_w_san} {move1_b_san}")
            print(f"2. {move2_w_san}")

            print("\nNumber of moves for White to mate: 2")
        else:
            print("\nSomething went wrong, the sequence did not result in a checkmate.")
            print("Current board:")
            print(board)

    except (ValueError, chess.InvalidMoveError, chess.IllegalMoveError) as e:
        print(f"\nAn error occurred while executing the moves: {e}")
        print("Please double-check the position and the move sequence.")

if __name__ == '__main__':
    find_mate_in_two()
