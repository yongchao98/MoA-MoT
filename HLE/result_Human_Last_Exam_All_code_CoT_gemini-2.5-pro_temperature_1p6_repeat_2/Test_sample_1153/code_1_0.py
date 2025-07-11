import chess
import chess.engine
import sys

def solve_chess_puzzle():
    """
    Analyzes a chess position with Stockfish to find the best move among a list of candidates.
    """
    # Define the chess position using Forsyth-Edwards Notation (FEN)
    fen = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"
    board = chess.Board(fen)

    # Define the candidate moves provided in the answer choices
    candidate_moves = {
        "A": {"label": "a8=Q", "uci": "a7a8q"},
        "B": {"label": "Nc5", "uci": "d3c5"},
        "C": {"label": "Kd4", "uci": "e3d4"},
        "D": {"label": "Kf2", "uci": "e3f2"},
        "E": {"label": "Nf4", "uci": "d3f4"},
        "F": {"label": "b4", "uci": "b3b4"}
    }
    
    # Locate the Stockfish engine.
    # Common paths are '/usr/games/stockfish' on Linux or just 'stockfish' if in PATH.
    # On Windows, you might need to provide the full path to 'stockfish.exe'.
    try:
        engine = chess.engine.SimpleEngine.popen_uci("stockfish")
    except FileNotFoundError:
        print("Error: Stockfish engine not found.", file=sys.stderr)
        print("Please install Stockfish and ensure it's in your system's PATH.", file=sys.stderr)
        return

    print(f"Analyzing position: {board.fen()}")
    print("White to move.")
    print("-" * 20)

    best_move_choice = None
    # Use a very small number for the initial best score
    best_score_cp = -float('inf') 
    results = []

    for choice, move_info in candidate_moves.items():
        move_label = move_info["label"]
        uci_move = move_info["uci"]
        
        move = chess.Move.from_uci(uci_move)
        
        # Make the move on a copy of the board to analyze the resulting position
        temp_board = board.copy()
        temp_board.push(move)
        
        # Analyze the position for a short duration
        analysis = engine.analyse(temp_board, chess.engine.Limit(time=0.5))
        
        # Get the score from White's perspective
        score = analysis["score"].white()
        
        # Format the score for display
        if score.is_mate():
            score_cp = 10000 * (1 if score.mate() > 0 else -1)
            score_str = f"Mate in {abs(score.mate())}"
        else:
            score_cp = score.score()
            # The evaluation is in centipawns, convert to pawn units
            score_str = f"{score_cp / 100.0:+.2f}"
            
        results.append((choice, move_label, score_str, score_cp))

    # Sort results from best to worst for white
    results.sort(key=lambda x: x[3], reverse=True)

    print("Evaluation of candidate moves (from best to worst for White):")
    for choice, label, score_str, score_cp in results:
        print(f"({choice}) {label}: {score_str}")
        
    best_move_choice, best_label, best_score_str, _ = results[0]
    
    print("-" * 20)
    print(f"The best move is {best_label}, with an evaluation of {best_score_str}.")
    
    # Close the engine process
    engine.quit()

    # Final answer in the required format
    print(f"<<<{best_move_choice}>>>")

if __name__ == "__main__":
    solve_chess_puzzle()