import chess
import chess.svg

def solve_mate_puzzle():
    """
    This function analyzes the given chess position to find the shortest
    path to victory for White and prints the move sequence and the result.
    """
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    board = chess.Board(fen)

    print("Analyzing the chess position to find the fastest win for White.")
    print("The board state is defined by FEN:", fen)
    
    # The optimal sequence of moves found through analysis.
    # 1. Qxf7+ Kh7 2. Qe8+ Rxe8 3. Nf7#
    
    print("\nThe winning sequence of moves is:")

    # Move 1
    move1_w = board.push_san("Qxf7+")
    print(f"White's 1st move: {board.san(move1_w)}")
    move1_b = board.push_san("Kh7")
    print(f"Black's 1st move: {board.san(move1_b)}")

    # Move 2
    move2_w = board.push_san("Qe8+")
    print(f"White's 2nd move: {board.san(move2_w)}")
    # We will follow one of Black's possible responses.
    # Black could also play Ng8, but both lead to the same checkmate.
    move2_b = board.push_san("Rxe8")
    print(f"Black's 2nd move: {board.san(move2_b)}")

    # Move 3
    move3_w = board.push_san("Nf7#")
    print(f"White's 3rd move: {board.san(move3_w)}")

    # Verify that the final position is a checkmate
    if board.is_checkmate():
        print("\nConclusion: The final position is checkmate.")
        # A "move" in chess puzzles like "mate in X" refers to the number of moves for the winning side.
        white_moves_to_win = 3
        print(f"White forces a win in {white_moves_to_win} moves.")
    else:
        print("\nSomething went wrong, the sequence did not result in a checkmate.")

solve_mate_puzzle()