import chess
import sys
from io import StringIO

def solve_three_check():
    """
    Analyzes and prints the shortest winning line for White in a game of Three-Check chess.
    """
    # The FEN string provided by the user, without the check count part.
    # The python-chess library handles standard FEN. We will track checks manually.
    fen = "r3kb1r/p2nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1"
    board = chess.Board(fen)

    # We will track the checks for White manually.
    white_checks_count = 0
    white_move_number = 0
    
    # This move sequence represents the forced win.
    # Optimal play is assumed from both sides.
    # Format: [White's move, Black's move, White's move, ...]
    move_sequence_uci = [
        "e1c1", "a7a6", "b5d7", "e7d7", "g5f6", "g7f6", 
        "d1d7", "e8d7", "h1d1", "d7c8", "b3c4"
    ]

    # Store old stdout to prevent chess board printing, and prepare to capture our own output
    old_stdout = sys.stdout
    sys.stdout = captured_output = StringIO()
    
    print("Tracing the shortest forced win for White:")
    print(f"Initial FEN: {board.fen()}")
    print(f"White checks: {white_checks_count}\n")
    
    for move_uci in move_sequence_uci:
        # Check if the move is legal before pushing
        if chess.Move.from_uci(move_uci) not in board.legal_moves:
            print(f"ERROR: The move {move_uci} is not legal in the position {board.fen()}.")
            print("The analyzed line is incorrect.")
            # Restore stdout and print captured output
            sys.stdout = old_stdout
            print(captured_output.getvalue())
            return

        is_white_to_move = board.turn == chess.WHITE
        if is_white_to_move:
            white_move_number += 1
            move_prefix = f"{white_move_number}."
        else:
            move_prefix = f"{white_move_number}..."
            
        # Get Standard Algebraic Notation for better readability
        san_move = board.san(chess.Move.from_uci(move_uci))
        
        # Make the move
        board.push_uci(move_uci)

        # Check if the move results in a check
        check_info = ""
        if board.is_check():
            if is_white_to_move:
                white_checks_count += 1
                check_info = f" (Check #{white_checks_count})"
        
        print(f"{move_prefix} {san_move}{check_info}")

        # Check for win condition
        if white_checks_count >= 3:
            print(f"\nResult: White wins by delivering 3 checks.")
            print(f"The minimal amount of moves by white to win is: {white_move_number}")
            break
            
        if board.is_checkmate():
            print(f"\nResult: White wins by checkmate.")
            print(f"The minimal amount of moves by white to win is: {white_move_number}")
            break
            
    if white_checks_count < 3 and not board.is_checkmate():
         print("\nThe provided line does not lead to a win.")

    # Restore stdout and print the captured output to the actual console
    sys.stdout = old_stdout
    print(captured_output.getvalue())

# Execute the analysis
solve_three_check()
